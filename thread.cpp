#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <time.h>


template <typename T>
class Matrix
{
private:
    int N, M;
    std::vector<std::vector<T> > data;
    
public:
    // Initialization from file
    Matrix(const std::string& path)
    {
        std::ifstream input(path);
        input >> N >> M;
        data.resize(N, std::vector<T>(M));
        for (int i = 0; i < N; ++i) 
        {
            for (int j = 0; j < M; ++j)
                input >> data[i][j];
        }
    }
    // Intialization with dimensions
    Matrix(int N, int M)
    {
        this->N = N;
        this->M = M;
        std::vector<T> row(M, 0);
        for (int i = 0; i < N; ++i)
            data.push_back(row);
    }

    // Copy initialization
    Matrix(const Matrix& other)
    {
        this->N = other.N;
        this->M = other.M;
        this->data = other.data;
    }

    // Print matrix
    void print() const
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < M; ++j)
                std::cout << data[i][j] << " ";
            std::cout << std::endl;
        }
    }

    // Write matrix to a file
    void write(const std::string& path) const
    {
        std::ofstream output(path);
        output << N << " " << M;
        for (int i = 0; i < N; ++i) 
        {
            output << std::endl;
            for (int j = 0; j < M; ++j)
                output << data[i][j] << " ";
        }
    }

    int get_rows_count() const {return N;}


    int get_cols_count() const {return M;}

    // Returns values from all cells
    std::vector<std::vector<T> > get_data() const {return data;}

    // Returns value of a single cell
    T get_item(int i, int j) const {return data[i][j];}

    // Changes value in single cell
    void set(int i, int j, T value)
    {
        data[i][j] = value;
    }

    // Transposes matrix
    void transpose() 
    {   
        if (N != M) throw std::invalid_argument("Only square matrices can be transposed.");
        T temp;
        for (int i = 0; i < N; ++i)
        {
            for (int j = i; j < N; ++j)
            {   
                temp = get_item(i, j);
                set(i, j, get_item(j, i));
                set(j, i, temp);
            }
        }
    }

    // Returns minor of a selected cell.
    double get_minor(int r, int c)
    {
        Matrix newm(N - 1, M - 1); 
        int t = 0;
        for (int i = 0; i < N; ++i) 
        {
            for (int j = 0; j < M; ++j)
            {
                if (i == r || j == c) continue;
                newm.set(t / (N - 1), t % (N - 1), get_item(i, j));
                ++t;                        
            }
        }
        return std::pow(-1, r + c) * get_det(newm);
    }

    // Returns zero size static matrix.
    static const Matrix& staticMatrix() {
        static Matrix<T> m(0, 0);
        return m;
    }

    // Returns determinant of this matrix
    double get_det(const Matrix& m = Matrix::staticMatrix()) 
    { 
        if (m.N == 0) return get_det(*this);
        if (m.N == 1)
            return m.get_item(0, 0);
        if (m.N == 2)
            return m.get_item(0, 0) * m.get_item(1, 1) - m.get_item(0, 1) * m.get_item(1, 0);
        else
        {
            double d = 0;
            for (int k = 0; k < m.N; ++k)
            {
                Matrix newm(m.N - 1, m.M - 1); 
                for (int i = 1; i < m.N; ++i) 
                {
                    int t = 0;
                    for (int j = 0; j < m.M; ++j)
                    {
                        if (j == k) continue;
                        newm.set(i - 1, t, m.get_item(i, j));
                        ++t;
                        
                    }
                }
                d += std::pow(-1, k + 2) * m.get_item(0, k) * get_det(newm);
            }
            return d;
        }
    }

    // Returns determinant of this matrix (THREADED)
    double get_det_threaded() 
    {
        if (N == 1) return get_item(0, 0);
        if (N == 2) return get_item(0, 0) * get_item(1, 1) - get_item(0, 1) * get_item(1, 0);
        else
        {
            double d = 0;
            std::vector<std::thread> threads;
            std::vector<double> determinants(N, 0);

            for (int k = 0; k < N; ++k)
            {
                threads.emplace_back([&, k]() 
                {
                    Matrix newm(N - 1, M - 1);
                    for (int i = 1; i < N; ++i)
                    {
                        int t = 0;
                        for (int j = 0; j < M; ++j)
                        {
                            if (j == k) continue;
                            newm.set(i - 1, t, get_item(i, j));
                            ++t;
                        }
                    }
                determinants[k] = std::pow(-1, k + 2) * get_item(0, k) * get_det(newm);
                });
        }

        for (auto& thread : threads)
        {
            thread.join();
        }

        for (double determinant : determinants)
        {
            d += determinant;
        }

        return d;
        }
    }


    // Returns product of substraction by other matrix
    Matrix<T> operator-(const Matrix& other) const
    {

        if (N != other.N or M != other.M)
        {
            throw std::invalid_argument("Matrices must have same sizes.");
        }
        
        Matrix NewMatrix(*this);
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < M; ++j)
                NewMatrix.data[i][j] -= other.data[i][j];
        } 

        return NewMatrix;
    }

    // Returns product of addition to other matrix
    Matrix<T> operator+(const Matrix& other) const
    {
        if (N != other.N or M != other.M)
        {
            throw std::invalid_argument("Matrices must have same sizes.");
        }
        
        Matrix NewMatrix(*this);
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < M; ++j)
                NewMatrix.data[i][j] += other.data[i][j];
        }
        return NewMatrix;
    }

    // Return product of multiplication by other matrix
    Matrix<T> operator*(const Matrix& other) const
    {
        if (N != other.M)
        {
            throw std::invalid_argument("Left matrix must have same columns count as left's rows count");
        }
        
        Matrix NewMatrix(N, other.M);
        int cell;
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < other.M; ++j)
            {
                cell = 0;
                for (int k = 0; k < M; ++k)
                {
                    cell += data[i][k] * other.data[k][j];
                }
                NewMatrix.set(i, j, cell);
            }   
        }
        return NewMatrix;
    }

    // Returns matrix multiplied by "a"
    Matrix<T> operator*(const int a) const
    {
        
        Matrix NewMatrix(*this);
        int cell;
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                NewMatrix.set(i, j, data[i][j] * a);
            }   
        }
        return NewMatrix;
    }

     // Returns matrix divided by integer
    Matrix<T> operator/(const int a) const
    {
        
        Matrix NewMatrix(*this);
        int cell;
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                NewMatrix.set(i, j, data[i][j] / a);
            }   
        }
        return NewMatrix;
    }

    // Comprasion with other matrix
    bool operator==(const Matrix& other) const 
    {
        if (N != other.N or M != other.M) {return false;}

        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < M; ++j)
                if (data[i][j] != other.data[i][j]) {return false;}
        }
        return true;
    }

    // Comprasion with scalar
    bool operator==(int scalar) const 
    {
        if (scalar == 1 || scalar == 0)
        {
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < M; ++j)
                {
                    switch (data[i][j])
                    {
                        case 0:
                            if (scalar && i == j) return false;
                            break;
                        case 1:
                            if (scalar && i != j) return false;
                            if (!scalar) return false;
                            break;
                        default:
                            return false;
                    }
                }
            }
            return true;
        }
        throw std::invalid_argument("Matrices can be compraised only with 0, 1 and other matrices.");
    }

    // Switching with other matrix
    void operator=(const Matrix& m)
    {
        this.N = m.get_rows_count();
        this.M = m.get_cols_count();
        this.data = m.get_data();
    }

    // Returns inverted matrix. (THREADED)
    Matrix<T> operator!()
    {
        double det = get_det_threaded();
        if (!det) throw std::runtime_error("Matrix cannot be reversed.");

        Matrix newm(N, M);
        transpose();

        std::vector<std::thread> threads;

        std::vector<std::vector<T> > minors(N, std::vector<T>(N, 0));

        for (int i = 0; i < N; ++i)
        {
            threads.emplace_back([&, i]() 
            {
                for (int j = 0; j < N; ++j)
                {
                    minors[i][j] = get_minor(i, j);
                }
            });
        }

        for (auto& thread : threads)
        {
            thread.join();
        }
        
        newm.data = minors;


        transpose();
        return newm / det;
    }
};
 

// Create identity matrix.
static Matrix<int> create_id_matrix(int a) 
{
    static Matrix<int> m(a, a);
    for (int i = 0; i < a; ++i) 
    {
        for (int j = 0; j < a; ++j) 
        {
            if (i == j) m.set(i, j, 1);
        }
    }
    return m;
}

// Create null matrix
static Matrix<int> create_null_matrix(int a) 
{
    static Matrix<int> m(a, a);
    return m;
}

