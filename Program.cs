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
            const string PATH = "txt\\2\\"; //путь, где лежат исходние данные 

            Node node = new Node();
            Element element = new Element();
            Area area = new Area();
            Boundary boundary = new Boundary();

            #region Чтение данных
            node.Read(PATH + "nodes.txt");
            area.Read(PATH + "area.txt");
            element.Read(PATH + "elements.txt", node, area);
            boundary.Read(PATH + "kryslov.txt", node, area);
            #endregion

            #region Вывод начальных данных (тест)
            //Console.WriteLine("********* Входные данные *******************");
            //Console.WriteLine("node");
            //node.Print();
            //Console.WriteLine();
            //Console.WriteLine("element");
            //element.Print();
            //Console.WriteLine();
            //Console.WriteLine("area");
            //area.Print();
            //Console.WriteLine();
            //material.Print();
            //Console.WriteLine();
            //Console.WriteLine("boundary");
            //boundary.Print();
            //Console.WriteLine("****************************");
            #endregion


            MKE(node, element, area, boundary);

            Console.ReadLine();
        }
        static Vector mat_mul(Vector ggl, Vector ggu, Vector ig, Vector jg, Vector di, Vector b)
        {
            Vector res = new Vector(); for (int i = 0; i < b.Items.Count; i++) res.Add();

            for (int i = 0; i < b.Items.Count(); i++)
            {
                res.Items[i].m = di.Items[i].m * b.Items[i].m;
                for (int j = ig.Items[i].a; j < ig.Items[i + 1].a; j++)
                {
                    res.Items[i].m += ggl.Items[j].m * b.Items[jg.Items[j].a].m;
                    res.Items[jg.Items[j].a].m += ggu.Items[j].m * b.Items[i].m;
                }
            }
            return res;
        }
        static Vector mat_mul_d(Vector di, Vector b)
        {
            Vector res = new Vector(); for (int i = 0; i < b.Items.Count; i++) res.Add();

            for (int i = 0; i < res.Items.Count; i++)
            {
                res.Items[i].m += di.Items[i].m * b.Items[i].m;
            }
            return res;
        }

        /// <summary>
        /// ЛОС
        /// </summary>
        static Vector loc(Vector ggl, Vector ggu, Vector ig, Vector jg, Vector di, Vector b)
        {
            double eps = 1e-20;
            int max_iter = 1000;

            for (int i = 0; i < di.Items.Count; i++)
            {
                di.Items[i].m = 1.0 / Math.Sqrt(di.Items[i].m);
            }
            for (int i = 0; i < ggl.Items.Count; i++)
            {
                ggl.Items[i].m = 0.0;
                ggu.Items[i].m = 0.0;
            }

            Vector x = new Vector(); for (int i = 0; i < b.Items.Count; i++) x.Add();

            Vector mat_mul_res = mat_mul(ggl, ggu, ig, jg, di, x);

            for (int i = 0; i < di.Items.Count; i++)
            {
                mat_mul_res.Items[i].m = b.Items[i].m - mat_mul_res.Items[i].m;
            }

            Vector r = new Vector(); for (int i = 0; i < mat_mul_res.Items.Count; i++) x.Add();
            r = mat_mul_d(di, mat_mul_res);

            double residual = r * r;

            double residual_next = residual + 1.0;
            double abs_residual_diff = Math.Abs(residual - residual_next);

            Vector z = new Vector(); for (int i = 0; i < r.Items.Count; i++) z.Add();
            z = mat_mul_d(di, r);

            mat_mul_res = mat_mul(ggl, ggu, ig, jg, di, z);

            Vector p = new Vector(); for (int i = 0; i < mat_mul_res.Items.Count; i++) p.Add();
            p = mat_mul_d(di, mat_mul_res);

            int k = 1;
            while (residual > eps && k < max_iter && abs_residual_diff > eps)
            {
                double pp = p * p;
                double alpha = (p * r) / pp;

                abs_residual_diff = Math.Abs(residual - residual_next);

                residual_next = residual;
                residual -= alpha * alpha * pp;

                for (int i = 0; i < di.Items.Count; i++)
                {
                    x.Items[i].m += alpha * z.Items[i].m;
                    r.Items[i].m -= alpha * p.Items[i].m;
                }

                Vector ur = new Vector(); for (int i = 0; i < r.Items.Count; i++) ur.Add();
                ur = mat_mul_d(di, r);
                mat_mul_res = mat_mul(ggl, ggu, ig, jg, di, ur);
                Vector dot_Rhs = mat_mul_d(di, mat_mul_res);

                double beta = -(p * dot_Rhs) / pp;

                for (int i = 0; i < di.Items.Count; i++)
                {
                    z.Items[i].m = ur.Items[i].m + beta * z.Items[i].m;
                    p.Items[i].m = dot_Rhs.Items[i].m + beta * p.Items[i].m;
                }

                k++;
            }

            return x;

        }

        static Matrix buildMatrix(Vector ggl, Vector ggu, Vector ig, Vector jg, Vector di)
        {
            Matrix res = new Matrix(); for (int i = 0; i < di.Items.Count; i++) res.Add();
            for (int i = 0; i < res.Items.Count; i++)
                for (int j = 0; j < res.Items.Count; j++) res.Items[i].r.Add();

            for (int i = 0; i < di.Items.Count; i++) res.Items[i].r.Items[i].m = di.Items[i].m;

            for (int i = 0; i < ig.Items.Count - 2; i++)
            {
                int k = ig.Items[i].a; //начальный элемент i-й строки в jg и ggl
                int n = ig.Items[i + 1].a - ig.Items[i].a; //количество элементов в i-й строке
                for (int j = 0; j < n; j++)
                {
                    res.Items[i].r.Items[jg.Items[k + j].a].m = ggl.Items[k + j].m;
                    res.Items[jg.Items[k + j].a].r.Items[i].m = ggu.Items[k + j].m;
                }
            }
            return res;
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

            #region Вычисление массивов ig, jg (портрет матрицы)
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
            ig.Add(ig.Items[ig.Items.Count - 1].a + 1);
            #endregion

            #region Нахождение матриц 
            Vector ggl = new Vector(); for (int i = 0; i < ig.Items[ig.Items.Count - 1].a - 1; i++) ggl.Add();
            Vector ggu = new Vector(); for (int i = 0; i < ig.Items[ig.Items.Count - 1].a - 1; i++) ggu.Add();
            Vector di = new Vector(); for (int i = 0; i < node.Items.Count; i++) di.Add();
            Vector b = new Vector(); for (int i = 0; i < node.Items.Count; i++) b.Add();

            foreach (Element item in element.Items)
            {
                item.createLocal();  // создание локальной матрицы и вектора (el.localMatrix и el.localVector)

                #region Глобальная матрица
                for (int j = 0; j < item.localMatrix.Items.Count; j++)
                {
                    di.Items[item.p[j]].m += item.localMatrix.Items[j].r.Items[j].m;
                    int i_beg = ig.Items[item.p[j]].a;
                    for (int k = 0; k < j; k++)
                    {
                        int i_end = ig.Items[item.p[j] + 1].a - 1;
                        while (jg.Items[i_beg].a != item.p[k])
                        {
                            int ind = (i_beg + i_end) / 2;
                            if (jg.Items[ind].a < item.p[k])
                            {
                                i_beg = ind + 1;
                            }
                            else
                            {
                                i_end = ind;
                            }
                        }
                        while (ggl.Items.Count <= i_beg) ggl.Add();
                        while (ggu.Items.Count <= i_beg) ggu.Add();
                        ggl.Items[i_beg].m += item.localMatrix.Items[j].r.Items[k].m;
                        ggu.Items[i_beg].m += item.localMatrix.Items[k].r.Items[j].m;
                        i_beg++;
                    }
                }

                //правая часть//
                for (int j = 0; j < 3; j++)
                {
                    b.Items[item.p[j]].m += item.localVector.Items[j].m;
                }
                #endregion
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

            foreach (Boundary item in boundary.Items)
            {
                item.setSecondCondition();
                item.setThirdCondition();

                for (int j = 0; j < item.localMatrix.Items.Count; j++)
                {
                    di.Items[item.p[j]].m += item.localMatrix.Items[j].r.Items[j].m;
                    int i_beg = ig.Items[item.p[j]].a;
                    for (int k = 0; k < j; k++)
                    {
                        int i_end = ig.Items[item.p[j] + 1].a - 1;
                        while (jg.Items[i_beg].a != item.p[k])
                        {
                            int ind = (i_beg + i_end) / 2;
                            if (jg.Items[ind].a < item.p[k])
                            {
                                i_beg = ind + 1;
                            }
                            else
                            {
                                i_end = ind;
                            }
                        }
                        while (ggl.Items.Count <= i_beg) ggl.Add();
                        while (ggu.Items.Count <= i_beg) ggu.Add();
                        ggl.Items[i_beg].m += item.localMatrix.Items[j].r.Items[k].m;
                        ggu.Items[i_beg].m += item.localMatrix.Items[k].r.Items[j].m;
                        i_beg++;
                    }
                }

                //правая часть//
                for (int j = 0; j < 2; j++)
                {
                    b.Items[item.p[j]].m += item.localVector.Items[j].m;
                }

            }


            Console.WriteLine("AFTER BOUND CONDITIONS (2,3)");

            Console.WriteLine("Global vector b:");
            b.Print();

            Console.WriteLine("Global matrix (ggl):");
            ggl.Print();

            Console.WriteLine("Global matrix (ggu):");
            ggu.Print();

            Console.WriteLine("DI:");
            di.Print();

            Matrix A = buildMatrix(ggl, ggu, ig, jg, di);
            Console.WriteLine("Матрица A:");
            A.Print();

            foreach (Boundary item in boundary.Items)
            {
                if (item.setFirstCondition())
                {
                    for (int k = 0; k < item.localMatrix.Items.Count; k++)
                    {
                        int i = item.p[k];
                        for (int j = 0; j < A.Items.Count; j++)
                        {
                            if (i == j)
                            {
                                A.Items[i].r.Items[j].m = 1;
                            }
                            else
                            {
                                A.Items[i].r.Items[j].m = 0;
                            }
                        }

                    }
                    //правая часть//
                    for (int j = 0; j < 2; j++)
                    {
                        b.Items[item.p[j]].m = item.localVector.Items[j].m;
                    }
                }
            }

            Console.WriteLine("Матрица A (AFTER BOUND CONDITIONS 1 )");
            A.Print();
            Vector gX = A / b;
            Console.WriteLine("Vector X:");
            gX.Print();

            //foreach (Boundary item in boundary.Items)
            //{
            //    if (item.setFirstCondition())
            //    {
            //        for (int j = 0; j < item.localMatrix.Items.Count; j++)
            //        {
            //            di.Items[item.p[j]].m = 1;
            //            int i_beg = ig.Items[item.p[j]].a;                        
            //            for (int k = 0; k < j; k++)
            //            {
            //                int i_end = ig.Items[item.p[j] + 1].a;
            //                for (int n= i_beg; n < i_end; n++)
            //                {
            //                    ggl.Items[n].m = 0;
            //                }
            //            }
            //        }

            //        //правая часть//
            //        for (int j = 0; j < 2; j++)
            //        {
            //            b.Items[item.p[j]].m = item.localVector.Items[j].m;
            //        }
            //    }
            //}

            Console.WriteLine("AFTER BOUND CONDITIONS 1 ");

            Console.WriteLine("Global vector b:");
            b.Print();

            //Console.WriteLine("Global matrix (ggl):");
            //ggl.Print();

            //Console.WriteLine("Global matrix (ggu):");
            //ggu.Print();

            //Console.WriteLine("DI:");
            //di.Print();

            //Vector gX = loc(ggl, ggu, ig, jg, di, b);
            //Console.WriteLine("Vector X:");
            //gX.Print();


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
                    Console.Write($"{ (String.Format("{0:f1}", Items[i].m))} ");
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
                        Console.Write($"{ (String.Format("{0:f1}", Items[i].r.Items[j].m))} ");
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
            static public Matrix operator *(double m, Matrix matr) //скалярное произведение матрицы на число
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

            static public Vector operator /(Matrix matr, Vector b) // A * x = b 
            {
                int n = matr.Items.Count;
                Vector x = new Vector(); for (int i = 0; i < b.Items.Count; i++) x.Add(b.Items[i].m);

                for (int i = 1; i < n; ++i)
                {
                    double sum = x.Items[i].m;
                    for (int j = 0; j < i; ++j)
                        sum -= matr.Items[i].r.Items[j].m * x.Items[j].m;
                    x.Items[i].m = sum;
                }
                x.Items[n - 1].m /= matr.Items[n - 1].r.Items[n - 1].m;
                for (int i = n - 2; i >= 0; --i)
                {
                    double sum = x.Items[i].m;
                    for (int j = i + 1; j < n; ++j)
                        sum -= matr.Items[i].r.Items[j].m * x.Items[j].m;
                    x.Items[i].m = sum / matr.Items[i].r.Items[i].m;
                }
                return x;
            }
        }

        /// <summary>
        /// Узлы (точки с координатами x,y)
        /// </summary>
        class Node
        {
            public double x, y;
            public List<Node> Items = new List<Node>();
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
                    Console.Write($"({Items[i].x}, {Items[i].y})\n");
                }
            }
        }
        /// <summary>
        /// Элементы (треугольники)
        /// </summary>
        class Element
        {
            public int[] p; //номера узлов
            public int a; // номер области
            public Node[] nodes; //узлы
            public Area area; // область в которой находится элемент
            public Matrix localMatrix = new Matrix(); //локальная мартрица элемента
            public Vector localVector = new Vector(); //локальная правая часть СЛАУ
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
            public void Print()
            {
                for (int i = 0; i < Items.Count; i++)
                {
                    Console.Write($"({Items[i].p[0]},{Items[i].p[1]},{Items[i].p[2]})-{Items[i].a} ");
                }
                Console.WriteLine();
            }
            public double detD() // detD (x1-x0)*(y2-y0) - (x2-x0)(y1-y0)
            {
                return (nodes[1].x - nodes[0].x) * (nodes[2].y - nodes[0].y) -
                    (nodes[2].x - nodes[0].x) * (nodes[1].y - nodes[0].y);
            }
            public double gamma(int k)
            {
                return 0; // nodes[k].y + nodes[k].x; //(формула x+y)
            }
            public double func(int k)
            {
                switch (this.a)
                {
                    case 0:
                        return -20;
                    case 1:
                        return 0;
                }
                return 2 * nodes[k].y; //какая-то функция для вычисления вектора правой части СЛАУ (например, 2*y)
            }
            public void createLocal() // Создание локальной матрици и вектора (localMatrix и localVector)
            {
                double detD = this.detD();
                double coef = Math.Abs(detD) / 60; //коэффициент для матрицы массы

                Matrix M = new Matrix(); //М.массы
                Matrix G = new Matrix(); //М.жескости

                for (int j = 0; j < 3; j++) //заполнение матрицы массы
                {
                    while (M.Items.Count <= j) M.Add();
                    for (int k = 0; k < 3; k++)
                    {
                        while (M.Items[j].r.Items.Count <= k) M.Items[j].r.Add();
                        double gamma = 0;

                        if (j == k)
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                gamma += (n == k ? 3 : 1) * this.gamma(n);
                            }
                        }
                        else
                        {
                            for (int n = 0; n < 3; n++)
                            {
                                gamma += this.gamma(n) / (n != j && n != k ? 2 : 1);
                            }
                        }
                        M.Items[j].r.Items[k].m = coef * gamma;
                    }
                }

                //коэффициенты альфа
                double[] alpha1 = new double[3];
                alpha1[0] = (this.nodes[1].y - this.nodes[2].y) / detD;   // (y1-y2) / detD
                alpha1[1] = (this.nodes[2].y - this.nodes[0].y) / detD; // (y2-y0) / detD
                alpha1[2] = (this.nodes[0].y - this.nodes[1].y) / detD; // (y0-y1) / detD

                double[] alpha2 = new double[3];
                alpha2[0] = (this.nodes[2].x - this.nodes[1].x) / detD; // (x2-x1) / detD
                alpha2[1] = (this.nodes[0].x - this.nodes[2].x) / detD; // (x0-x2) / detD
                alpha2[2] = (this.nodes[1].x - this.nodes[0].x) / detD; // (x1-x0) / detD


                for (int j = 0; j < 3; j++) //заполнение матрицы жесткости
                {
                    G.Add();
                    for (int k = 0; k < 3; k++)
                    {
                        G.Items[j].r.Add();
                        G.Items[j].r.Items[k].m = this.area.lymbda * Math.Abs(detD) * (alpha1[j] * alpha1[k] + alpha2[j] * alpha2[k]) / 2;
                    }
                }

                localMatrix = M + G;
                coef = detD / 24;
                localVector.Add(coef * (2 * this.func(0) + 1 * this.func(1) + 1 * this.func(2)));
                localVector.Add(coef * (1 * this.func(0) + 2 * this.func(1) + 1 * this.func(2)));
                localVector.Add(coef * (1 * this.func(0) + 1 * this.func(1) + 2 * this.func(2)));
                //localVector = M * localVector;


                //Console.WriteLine($"<<<<<<< {this.a} >>>>>>>>");
                //Console.WriteLine("Local M:");
                //M.Print();
                //Console.WriteLine("Local G:");
                //G.Print();

                //Console.WriteLine("Local A:");
                //localMatrix.Print();
                //Console.WriteLine("Local b:");
                //localVector.Print();
                //Console.WriteLine("--------------------------");

            }

        }
        /// <summary>
        /// Значения параметрлв для областей
        /// лямбда и бета
        /// </summary>
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
                    Console.Write($"lyambda={Items[i].lymbda} beta={Items[i].beta} ");
                    Console.WriteLine();
                }
            }
        }

        /// <summary>
        /// Краевые условия
        /// номер области, две границы ребра, типа краевого условия (1,2,3), номер условия (номер функции)
        /// </summary>
        class Boundary
        {
            public int a; // номер области
            public int[] p; //номера узлов границ ребра
            public int condType; //тип условия (1-е,2-е,3-е)
            public int condNum; //номер условия для выбора функции 
            public Node[] nodes; //узлы
            public Area area; // область в которой находится элемент
            public Matrix localMatrix = new Matrix(); //локальная мартрица элемента
            public Vector localVector = new Vector(); //локальная правая часть СЛАУ
            public List<Boundary> Items = new List<Boundary>();
            public void Read(string filename, Node node, Area _area)
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
                            p = new int[2] { int.Parse(words[1]), int.Parse(words[2]) },
                            condType = int.Parse(words[3]),
                            condNum = int.Parse(words[4]),
                            area = _area.Items[int.Parse(words[0])],
                            nodes = new Node[2] { node.Items[int.Parse(words[1])], node.Items[int.Parse(words[2])] }
                        });
                    }
                }
            }
            public void Print()
            {
                for (int i = 0; i < Items.Count; i++)
                {
                    Console.Write($"{Items[i].a},{Items[i].p[0]},{Items[i].p[1]},{Items[i].condType},{Items[i].condNum} ");
                    Console.WriteLine();
                }
            }
            public double mesG() //длина ребра
            {
                return Math.Sqrt(Math.Pow(nodes[0].x - nodes[1].x, 2.0) + Math.Pow(nodes[0].y - nodes[1].y, 2.0));
            }
            public double func1(int k)
            {
                return Math.Pow(nodes[k].y, 2);
            }
            public double func2(int k)
            {
                switch (this.condNum)
                {
                    case 0:
                        return 20;
                    case 1:
                        return 0;
                }
                return 0;
            }
            public double func3(int k)
            {
                return 20 * nodes[k].y - 27;
            }
            // 1-е краевое условие
            public bool setFirstCondition()
            {
                if (this.condType == 1)
                {
                    localVector.Items[0].m += this.func1(0);
                    localVector.Items[1].m += this.func1(1);

                    localMatrix.Items[0].r.Items[0].m = 1; localMatrix.Items[0].r.Items[1].m = 0;
                    localMatrix.Items[1].r.Items[0].m = 0; localMatrix.Items[1].r.Items[1].m = 1;

                }
                return (this.condType == 1);
            }
            // 2-е краевое условие
            public void setSecondCondition()
            {
                if (this.condType == 2)
                {
                    double coef = this.mesG() / 6;

                    localVector.Items[0].m += coef * (2 * this.func2(0) + 1 * this.func2(1));
                    localVector.Items[1].m += coef * (1 * this.func2(0) + 2 * this.func2(1));

                }
            }
            // 3-е краевое условие
            public void setThirdCondition()
            {
                if (this.condType == 3)
                {
                    double coef = this.area.beta * this.mesG() / 6;
                    localVector.Items[0].m += coef * (2 * this.func3(0) + 1 * this.func3(1));
                    localVector.Items[1].m += coef * (1 * this.func3(0) + 2 * this.func3(1));

                    localMatrix.Items[0].r.Items[0].m += (coef * 2); localMatrix.Items[0].r.Items[1].m += (coef * 1);
                    localMatrix.Items[1].r.Items[0].m += (coef * 1); localMatrix.Items[1].r.Items[1].m += (coef * 2);
                }
            }

            public Boundary()
            {
                while (localVector.Items.Count <= 1) localVector.Add();
                while (localMatrix.Items.Count <= 1) localMatrix.Add();
                while (localMatrix.Items[0].r.Items.Count <= 1) localMatrix.Items[0].r.Add();
                while (localMatrix.Items[1].r.Items.Count <= 1) localMatrix.Items[1].r.Add();
            }
        }
    }
}

