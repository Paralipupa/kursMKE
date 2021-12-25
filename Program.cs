using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace appMKE
{

    internal class Program
    {
        static void Main(string[] args)
        {
            const string PATH = "txt\\0\\"; //путь, где лежат исходние данные 

            Node node = new Node();
            Element element = new Element();
            Area area = new Area();
            Material material = new Material();
            Boundary boundary = new Boundary();

            #region Чтение данных
            node.Read(PATH + "nodes.txt");
            area.Read(PATH + "area.txt");
            element.Read(PATH + "elements.txt", node, area);
            material.Read(PATH + "mat.txt");
            boundary.Read(PATH + "kryslov.txt");
            #endregion

            #region Вывод начальных данных (тест)
            //Console.WriteLine("node");
            //node.Print();
            //Console.WriteLine();
            //Console.WriteLine("element");
            //element.Print();
            //Console.WriteLine();
            //Console.WriteLine("area");
            //area.Print();
            //Console.WriteLine();
            //Console.WriteLine("mat");
            //material.Print();
            //Console.WriteLine();
            //Console.WriteLine("boundary");
            //boundary.Print();
            #endregion


            MKE(node, element, area, boundary);

            Console.ReadLine();
        }
                
        static void MKE(Node node, Element element, Area area, Boundary boundary)
        {
            #region Построение списка связаности
            Dictionary<int, int>[] listNodes = new Dictionary<int, int>[node.Items.Count];
            for (int i = 0; i < listNodes.Length; i++)
            {
                listNodes[i] = new Dictionary<int, int>();
            }
            for (int i = 0; i < element.Items.Count; i++)
            {
                int[] inodes = { element.Items[i].p[0], element.Items[i].p[1], element.Items[i].p[2] };
                for (int j = 0; j < 2; j++)
                {
                    for (int k = j + 1; k < 3; k++)
                    {
                        if (inodes[j] > inodes[k])
                        {
                            int p = inodes[j];
                            inodes[j] = inodes[k];
                            inodes[k] = p;
                        }
                    }
                }

                if (listNodes[inodes[2]].ContainsKey(inodes[1]) == false)
                {
                    listNodes[inodes[2]].Add(inodes[1], inodes[1]);
                }
                if (listNodes[inodes[2]].ContainsKey(inodes[0]) == false)
                {
                    listNodes[inodes[2]].Add(inodes[0], inodes[0]);
                }
                if (listNodes[inodes[1]].ContainsKey(inodes[0]) == false)
                {
                    listNodes[inodes[1]].Add(inodes[0], inodes[0]);
                }
            }
            for (int i = 0; i < listNodes.Length; i++)
            {
                listNodes[i] = listNodes[i].OrderBy(x => x.Key).ToDictionary(x => x.Key, x => x.Value);
            }
            #endregion

            #region Вычисление массивов ig, jg
            Vector ig = new Vector();
            Vector jg = new Vector();

            ig.Add();
            ig.Add();

            for (int i = 0; i < listNodes.Length; i++)
            {
                foreach (KeyValuePair<int, int> keyVal in listNodes[i])
                {
                    jg.Add(keyVal.Key);
                }
            }

            for (int i = 2; i < node.Items.Count + 1; i++)
            {
                ig.Add(ig.Items[i - 1].a + listNodes[i - 1].Count);
            }
            #endregion

            #region Нахождение матриц 
            Matrix M = new Matrix(); //М.массы
            Matrix G = new Matrix(); //М.жескости
            Vector ggl = new Vector();
            Vector ggu = new Vector();
            Vector di = new Vector(); for (int i = 0; i < node.Items.Count; i++) di.Add();
            Vector b = new Vector(); for (int i = 0; i < node.Items.Count; i++) b.Add();

            foreach (Element el in element.Items)
            {
                double detD = (node.Items[el.p[1]].x - node.Items[el.p[0]].x) * (node.Items[el.p[2]].y - node.Items[el.p[0]].y) -
                    (node.Items[el.p[2]].x - node.Items[el.p[0]].x) * (node.Items[el.p[2]].y - node.Items[el.p[1]].y);
                double coef = Math.Abs(detD) / 60; //коэффициент для матрицы массы
                for (int j = 0; j < 3; j++) //заполнение матрицы массы
                {
                    M.Add();
                    for (int k = 0; k < 3; k++)
                    {
                        M.Items[j].r.Add();
                        double gamma = 0;

                        if (j == k)
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                gamma += (n == k ? 3 : 1) * el.gamma(n);
                            }
                        }
                        else
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                gamma += el.gamma(n) / (n != j && n != k ? 2 : 1);
                            }
                        }
                        M.Items[j].r.Items[k].m = coef * gamma;
                    }
                }

                Vector local_b = new Vector();

                local_b.Add(1.0); local_b.Add(1.0); local_b.Add(2.0); //здесь нужны данные вектора b

                local_b = M * local_b;

                for (int j = 0; j < 3; j++)
                {
                    b.Items[el.p[j]].m += local_b.Items[j].m;
                }

                //коэффициенты альфа
                double[] alpha1 = new double[3];
                alpha1[0] = (node.Items[el.p[1]].y - node.Items[el.p[2]].y) / detD;
                alpha1[1] = (node.Items[el.p[2]].y - node.Items[el.p[0]].y) / detD;
                alpha1[2] = (node.Items[el.p[0]].y - node.Items[el.p[1]].y) / detD;

                double[] alpha2 = new double[3];
                alpha2[0] = (node.Items[el.p[2]].x - node.Items[el.p[1]].x) / detD;
                alpha2[1] = (node.Items[el.p[0]].x - node.Items[el.p[2]].x) / detD;
                alpha2[2] = (node.Items[el.p[1]].x - node.Items[el.p[0]].x) / detD;

                for (int j = 0; j < 3; j++) //заполнение матрицы жесткости
                {
                    G.Add();
                    for (int k = 0; k < 3; k++)
                    {
                        G.Items[j].r.Add();
                        G.Items[j].r.Items[k].m = el.area.lymbda * Math.Abs(detD) * (alpha1[j] * alpha1[k] + alpha2[j] * alpha2[k]) / 2;
                    }
                }

                Console.WriteLine("Local M:");
                M.Print();
                Console.WriteLine("Local G:");
                G.Print();

                Matrix local_A = new Matrix();
                local_A = M + G; //локальная матрица А
                Console.WriteLine("Local A:");
                local_A.Print();

                #region Глобальная матрица

                //for (int j = 0; j < local_A.Items.Count; j++)
                //{
                //    di.Items[el.p[j]].m += local_A.Items[j].r.Items[j].m;
                //    int i_beg = ig.Items[el.p[j]].a;
                //    for (int k = 0; k < j; k++)
                //    {
                //        int i_end = ig[element.joints[j] + 1] - 1;
                //        while (jg[i_beg] != element.joints[k])
                //        {
                //            int ind = (i_beg + i_end) / 2;
                //            if (jg[ind] < element.joints[k])
                //            {
                //                i_beg = ind + 1;
                //            }
                //            else
                //            {
                //                i_end = ind;
                //            }
                //        }
                //        ggl[i_beg] += local_A[j][k];
                //        ggu[i_beg] += local_A[k][j];
                //        i_beg++;
                //    }
                //}

                #endregion


                /*      
                 *          add_to_global(ggl, ggu, di, ig, jg, local_A, element); //Глобальная матрица
                */
            }
            #endregion

            #region Вывод данных

            Console.WriteLine("BEFORE BOUND CONDITIONS");

            Console.WriteLine("Global vector b:");
            b.Print();

            Console.WriteLine("Global matrix (ggl):");
            ggl.Print();

            Console.WriteLine("Global matrix (ggu):");
            ggu.Print();

            Console.WriteLine("DI:");
            di.Print();

            //boundary_conditions(nodes, boundaries, areas, di, ggl, ggu, ig, jg, b);

            Console.WriteLine("AFTER BOUND CONDITIONS");

            Console.WriteLine("Global vector b:");
            b.Print();

            Console.WriteLine("Global matrix (ggl):");
            ggl.Print();

            Console.WriteLine("Global matrix (ggu):");
            ggu.Print();

            Console.WriteLine("DI:");
            di.Print();

            Vector gX = new Vector();
            //gX = loc(ggl, ggu, ig, jg, di, b);
            Console.WriteLine("Vector X:");
            gX.Print();

            #endregion
        }


        class Vector
        {
            public int a;
            public double m;
            public List<Vector> Items = new List<Vector>();
            public void Add()
            {
                Items.Add(new Vector { a = 0, m = 0 });
            }
            public void Add(int _a)
            {
                Items.Add(new Vector { a = _a, m = 0 });
            }
            public void Add(double _x)
            {
                Items.Add(new Vector { a = 0, m = _x });
            }
            public void Print()
            {
                for (int i = 0; i < Items.Count; i++)
                {
                    Console.Write($"{Items[i].m} ");
                }
                Console.WriteLine();
            }
            static public double operator *(Vector V1, Vector V2) //скалярное произведение векторов
            {
                double res = 0.0;
                for (int i = 0; i < V1.Items.Count; i++)
                {
                    res += V1.Items[i].m * V2.Items[i].m;
                }
                return res;
            }
            static public double Norm(Vector x)
            {
                return Math.Sqrt(x * x);
            }
        }

        class Matrix
        {
            public Vector r = new Vector();
            public List<Matrix> Items = new List<Matrix>();
            public void Add()
            {
                Items.Add(new Matrix());
            }
            public void Print()
            {
                for (int i = 0; i < Items.Count; i++)
                {
                    for (int j = 0; j < Items[i].r.Items.Count; j++)
                    {
                        Console.Write($"{Items[i].r.Items[j].m} ");
                    }
                    Console.WriteLine();
                }
            }
            static public Vector operator *(Matrix matr, Vector vect) //скалярное произведение векторова и матрицы
            {
                Vector res = new Vector();
                for (int i = 0; i < vect.Items.Count; i++)
                {
                    res.Add();
                    for (int j = 0; j < vect.Items.Count; j++)
                    {
                        res.Items[i].m += matr.Items[i].r.Items[j].m * vect.Items[j].m;
                    }
                }
                return res;
            }
            static public Matrix operator *(double m, Matrix matr) //скалярное произведение векторова и матрицы
            {
                Matrix res = new Matrix();
                for (int i = 0; i < matr.Items.Count; i++)
                {
                    res.Add();
                    for (int j = 0; j < matr.Items[i].r.Items.Count; j++)
                    {
                        res.Items[i].r.Add();
                        res.Items[i].r.Items[j].m += matr.Items[i].r.Items[j].m * m;
                    }
                }
                return res;
            }
            static public Matrix operator +(Matrix matrA, Matrix matrB) //сумма матриц
            {
                Matrix res = new Matrix();
                for (int i = 0; i < matrA.Items.Count; i++)
                {
                    res.Add();
                    for (int j = 0; j < matrA.Items[i].r.Items.Count; j++)
                    {
                        res.Items[i].r.Add();
                        res.Items[i].r.Items[j].m += (matrA.Items[i].r.Items[j].m + matrB.Items[i].r.Items[j].m);
                    }
                }
                return res;
            }
        }
        /// <summary>
        /// Узлы (точки с координатами x,y)
        /// </summary>
        class Node
        {
            public double x, y;
            public List<Node> Items = new List<Node>();
            public double det(Node node1, Node node2, Node node3)
            {
                return ((node2.x - node1.x) * (node3.y - node1.y) -
                    (node3.x - node1.x) * (node2.y - node1.y));
            }
            public void Read(string filename)
            {
                string[] lines = File.ReadAllLines(filename);
                foreach (string s in lines)
                {
                    string[] words = s.Split(new char[] { ' ' });
                    if (words.Length == 2)
                    {
                        Items.Add(new Node
                        {
                            x = double.Parse(words[0]),
                            y = double.Parse(words[1])
                        });
                    }
                }
            }
            public void Print()
            {
                for (int i = 0; i < Items.Count; i++)
                {
                    Console.Write($"{Items[i].x} {Items[i].y}\n");
                }
            }
        }
        /// <summary>
        /// Элементы (треугольники)
        /// </summary>
        class Element
        {
            public int[] p; //номера узлов
            public Node[] nodes; //узлы
            public int a; // номер области
            public Area area; // область в которой находится элемент
            public List<Element> Items = new List<Element>();
            public void Read(string filename, Node node, Area _area)
            {
                string[] lines = File.ReadAllLines(filename);
                foreach (string s in lines)
                {
                    string[] words = s.Split(new char[] { ' ' });
                    if (words.Length == 4)
                    {
                        Items.Add(new Element
                        {
                            p = new int[3] { int.Parse(words[0]), int.Parse(words[1]), int.Parse(words[2]) },
                            a = int.Parse(words[3]),
                            nodes = new Node[3] { node.Items[int.Parse(words[0])], node.Items[int.Parse(words[1])], node.Items[int.Parse(words[2])] },
                            area = _area.Items[int.Parse(words[3])]
                        });
                    }
                }
            }
            public double gamma(int k)
            {
                return nodes[k].y + nodes[k].x; //(формула x+y)
            }

            public void Print()
            {
                for (int i = 0; i < Items.Count; i++)
                {
                    Console.Write($"({Items[i].p[0]},{Items[i].p[1]},{Items[i].p[2]}) ");
                }
                Console.WriteLine();
            }
        }

        class boundary
        {
            public int area, node0, node1, condType, condNum;

        }

        class Area
        {
            public double lymbda, beta;
            public List<Area> Items = new List<Area>();
            public void Read(string filename)
            {
                string[] lines = File.ReadAllLines(filename);
                foreach (string s in lines)
                {
                    string[] words = s.Split(new char[] { ' ' });
                    if (words.Length == 2)
                    {
                        Items.Add(new Area
                        {
                            lymbda = int.Parse(words[0]),
                            beta = int.Parse(words[1])
                        });
                    }
                }
            }
            public void Print()
            {
                for (int i = 0; i < Items.Count; i++)
                {
                    Console.Write($"beta={Items[i].lymbda} gamma={Items[i].beta} ");
                }
                Console.WriteLine();
            }
        }

        class Material
        {
            List<List<double>> param = new List<List<double>>();
            List<double> fa = new List<double>();
            public void Read(string filename)
            {
                string[] lines = File.ReadAllLines(filename);
                foreach (string s in lines)
                {
                    string[] words = s.Split(new char[] { ' ' });
                    if (words.Length == 2)
                    {
                        param.Add(new List<double> { double.Parse(words[0]), double.Parse(words[1]) });
                    }
                    else if (words.Length == 1)
                    {
                        fa.Add(double.Parse(words[0]));
                    }
                }
            }
            public void Print()
            {
                for (int i = 0; i < param.Count; i++)
                {
                    Console.Write($"{param[i][0]} {param[i][1]} ");
                    Console.WriteLine();
                }
                for (int i = 0; i < fa.Count; i++)
                {
                    Console.Write($"{fa[i]} ");
                }
                Console.WriteLine();
            }
        }
        class Boundary
        {
            public int a, p0, p1, condType, condNum;
            List<Boundary> Items = new List<Boundary>();
            public void Read(string filename)
            {
                string[] lines = File.ReadAllLines(filename);
                foreach (string s in lines)
                {
                    string[] words = s.Split(new char[] { ' ' });
                    if (words.Length == 5)
                    {
                        Items.Add(new Boundary
                        {
                            a = int.Parse(words[0]),
                            p0 = int.Parse(words[1]),
                            p1 = int.Parse(words[2]),
                            condType = int.Parse(words[3]),
                            condNum = int.Parse(words[4])
                        });
                    }
                }
            }
            public void Print()
            {
                for (int i = 0; i < Items.Count; i++)
                {
                    Console.Write($"{Items[i].a},{Items[i].p0},{Items[i].p1},{Items[i].condType},{Items[i].condNum} ");
                    Console.WriteLine();
                }
            }

        }
    }
}

