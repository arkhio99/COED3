using System;
using System.Data.Common;
using System.Text;

namespace Matrix
{
    public class Matrix {
        private double[,] _val;
        private double _epsilon;
        public double this[Int32 i, Int32 j]
        {
            get {
                return _val[i, j];
            }
            set {
                _val[i, j] = value;
            }
        }
        public int rows {
            get {
                return _val.GetLength(0);
            }
        }
        public int cols {
            get {
                return _val.GetLength(1);
            }
        }
        public double Epsilon
        {
            get
            {
                return _epsilon;
            }
            set
            {
                _epsilon = value;
            }
        }
        public Matrix GetRow(int r)
        {
            double[,] row = new double[1, cols];
            for (int i = 0; i < cols; i++)
                row[0,i] = _val[r, i];
            return new Matrix(row, _epsilon);
        }

        public Matrix GetColumn(int c)
        {
            double[,] col = new double[rows, 1];
            for (int i = 0; i < cols; i++)
                col[i, 0] = _val[i, c];
            return new Matrix(col, _epsilon);
        }

        public void DeleteCol(int c)
        {
            double[,] newM = new double[rows, cols - 1];
            for (int i = 0 ;i<rows;i++)
            {
                int offset = 0;
                for(int j=0; j<cols;j++)
                {
                    if (j == c)
                        offset = -1;
                    else
                        newM[i, j + offset] = _val[i, j];
                }
            }
            _val = newM;
        }

        public void AddCol()
        {
            double[,] newM = new double[rows, cols + 1];
            for(int i=0;i<rows;i++)
                for(int j=0;j<cols;j++)
                {
                    newM[i, j] = _val[i, j];
                }
            for(int i=0;i<rows;i++)
            {
                newM[i, cols] = 0;
            }
            _val = newM;
        }
        public Matrix(double[,] ar, double eps=0.001)
        {
            _val=new double[ar.GetLength(0),ar.GetLength(1)];
            for (int i = 0; i < ar.GetLength(0); i++)
            {
                int n = ar.GetLength(1);
                for (int j = 0; j < n; j++)
                {
                    _val[i, j] = ar[i, j];
                }
            }
            _epsilon=eps;
        }
        public Matrix(int _r, int _c, double _eps=0.001)
        {
            _val = new double[_r, _c];
            _epsilon = _eps;
        }
        public Matrix(Matrix m)
        {
            int r = m.rows;
            int c = m.cols;
            _val = new double[r, c];
            for(int i=0;i<r;i++)
                for(int j=0;j<c;j++)
                {
                    _val[i, j] = m[i, j];
                }
            _epsilon = m.Epsilon;
        }
        public Matrix Transp()
        {
            double[,] res=new double[cols,rows];
            for(int i=0;i<rows;i++)
                for(int j=0;j<cols;j++)
                    res[j,i]=_val[i,j];
            return new Matrix(res,Epsilon);
        }

        private void SwapStrs(double[,] ar,int i, int j)
        {
            for(int k=0;k<cols;k++)
            {
                ar[i,k]+=ar[j,k];
                ar[j,k]=ar[i,k]-ar[j,k];
                ar[i,k]-=ar[j,k];
            }
        }
        public Matrix GetDiagMatr(out double mnoj)
        {
            mnoj=1;
            double[,] res=new double[rows,cols];
            for(int i=0;i<rows;i++)
                for(int j=0;j<cols;j++)
                    res[i,j]=_val[i,j];
            for(int i=0;i<rows;i++)
            {
                if(Math.Abs(res[i,i])<_epsilon)
                {
                    for(int j=0;j<rows;j++)
                    if(Math.Abs(res[j,i])>_epsilon)
                    {
                        SwapStrs(res,i,j);
                        mnoj*=-1;
                    }
                }
                for(int j=0;j<rows;j++)
                {
                    if(j==i || res[i,i]==0)
                    continue;
                    double temp=res[j,i];
                    for(int k=0;k<cols;k++)
                    {
                        res[j,k]-=res[i,k]/res[i,i]*temp;
                    }                    
                }                
            }
            return new Matrix(res,_epsilon);
        }
        public double GetDet()
        {
            double res=1;
            Matrix diag=GetDiagMatr(out double mnoj);
            for(int i=0;i<rows;i++)
            res*=diag[i,i];
            return res*mnoj;            
        }
        public Matrix GetRevMatr()
        {
            if (Math.Abs(GetDet()) < _epsilon)
                throw new Exception("Determinant equals zero");
            if (rows != cols)
                throw new Exception("Not square matrix");
            double[,] temp_ar = new double[rows, cols];
            double[,] obr = new double[rows, cols];
            for(int i=0;i<rows;i++)
                for(int j=0;j<cols;j++)
                {
                    temp_ar[i, j] = _val[i, j];
                    obr[i, j] = i == j ? 1 : 0;
                }
            for (int i = 0; i < rows; i++)
            {
                if (Math.Abs(temp_ar[i, i]) < _epsilon)
                {
                    for (int j = 0; j < rows; j++)
                        if (Math.Abs(temp_ar[j, i]) > _epsilon)
                        {
                            SwapStrs(temp_ar, i, j);
                            SwapStrs(obr, i, j);
                        }
                }
                for (int j = 0; j < rows; j++)
                {
                    if (j == i || temp_ar[i, i] == 0)
                        continue;
                    double temp = temp_ar[j, i];
                    for (int k = 0; k < cols; k++)
                    {
                        temp_ar[j, k] -= temp_ar[i, k] / temp_ar[i, i] * temp;
                        obr[j, k] -= obr[i, k] / temp_ar[i, i] * temp;
                    }

                }
            }
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    obr[i, j] /= temp_ar[i, i];
            return new Matrix(obr,_epsilon);
        }

        public override string ToString()
        {
            StringBuilder s=new StringBuilder();
            for(int i=0;i<rows;i++)
            {
                for(int j=0;j<cols;j++)
                {
                    s.Append(Math.Round(_val[i,j],3).ToString());
                    s.Append("\t");
                }
                s.Append("\n");
            }
            return s.ToString();
        }

        static public Matrix operator +(Matrix a, Matrix b)
        {
            if (a.rows != b.rows && b.cols != a.cols)
                throw new Exception("Sizes of matrixes not equal");
            double[,] newM =new double[a.rows, a.cols];
            for (int r = 0; r < a.rows; r++)
                for (int c = 0; c < a.cols; c++)
                    newM[r, c] = a[r, c] + b[r, c];
            return new Matrix(newM, a.Epsilon);
        }
        static public Matrix operator -(Matrix a, Matrix b)
        {
            if (a.rows != b.rows && b.cols != a.cols)
                throw new Exception("Sizes of matrixes not equal");
            double[,] newM = new double[a.rows, a.cols];
            for (int r = 0; r < a.rows; r++)
                for (int c = 0; c < a.cols; c++)
                    newM[r, c] = a[r, c] - b[r, c];
            return new Matrix(newM, a.Epsilon);
        }
        static public Matrix operator *(Matrix a, Matrix b)
        {
            if (a.cols != b.rows)
                throw new Exception("Number of columns of left matrix is not equal to number of rows of right matrix");
            double[,] newM = new double[a.rows, b.cols];
            for(int i=0;i<a.rows;i++)
                for(int j=0;j<b.cols;j++)
                {
                    newM[i, j] = 0;
                    for (int k = 0; k < a.cols; k++)
                        newM[i, j] += a[i, k] * b[k, j];
                }
            return new Matrix(newM, a.Epsilon);
        }
    }
}