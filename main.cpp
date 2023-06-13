#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>


template <typename T>
class Matrix
{
private:
    int N, M;
    std::vector<std::vector<T> > data;
    
public:
    // Initialization from file
    Matrix(const std::string& );

    // Intialization with dimensions
    Matrix(int, int);

    // Copy initialization
    Matrix(const Matrix&);

    // Writes matrix to a file
    void write(const std::string&) const;

    // Returns rows count
    int get_rows() const;

    // Returns columns count
    int get_cols() const;

    // Returns values from all cells
    std::vector<std::vector<T> > get_data() const;

    // Returns value of a single cell
    T get_item(int, int) const;

    // Changes value in single cell
    void set(int, int, T);
    
    // Transposes matrix
    void transpose();

    // Returns minor of a selected cell
    double get_minor(int, int);

    // Returns zero size static matrix
    static const Matrix& staticMatrix()
    {
    static Matrix<T> m(0, 0);
    return m;
    }

    // Returns determinant of this matrix
    double get_det(const Matrix& m=Matrix::staticMatrix());

    // Returns product of substraction by other matrix
    Matrix<T> operator-(const Matrix&) const;

    // Returns product of addition to other matrix
    Matrix<T> operator+(const Matrix&) const;

    // Return product of multiplication by other matrix
    Matrix<T> operator*(const Matrix&) const;

    // Returns matrix multiplied by integer
    Matrix<T> operator*(const int) const;

    // Returns matrix divided by integer
    Matrix<T> operator/(const int) const;

    // Returns result of comprasion with other matrix
    bool operator==(const Matrix&) const;

    // Returns result of comprasion with scalar
    bool operator==(int) const;

    // Switches with other matrix
    void operator=(const Matrix&);
    
    // Returns inverted matrix
    Matrix<T> operator!();
};


// Initialization from file
template <typename T>
Matrix<T>::Matrix(const std::string& path)
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
template <typename T>
Matrix<T>::Matrix(int N, int M)
{
    this->N = N;
    this->M = M;
    std::vector<T> row(M, 0);
    for (int i = 0; i < N; ++i)
        data.push_back(row);
}

// Copy initialization
template <typename T>
Matrix<T>::Matrix(const Matrix& other)
{
    this->N = other.N;
    this->M = other.M;
    this->data = other.data;
}

// Write matrix to a file
template <typename T>
void Matrix<T>::write(const std::string& path) const
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

template <typename T>
int Matrix<T>::get_rows() const {return N;}

template <typename T>
int Matrix<T>::get_cols() const {return M;}

// Returns values from all cells
template <typename T>
std::vector<std::vector<T> > Matrix<T>::get_data() const {return data;}

// Returns value of a single cell
template <typename T>
T Matrix<T>::get_item(int i, int j) const {return data[i][j];}

// Changes value in single cell
template <typename T>
void Matrix<T>::set(int i, int j, T value)
{
    data[i][j] = value;
}

// Transposes matrix
template <typename T>
void Matrix<T>::transpose() 
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
template <typename T>
double Matrix<T>::get_minor(int r, int c)
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

// Returns determinant of this matrix
template <typename T>
double Matrix<T>::get_det(const Matrix& m)
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

// Returns product of substraction by other matrix
template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix& other) const
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
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& other) const
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
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& other) const
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
template <typename T>
Matrix<T> Matrix<T>::operator*(const int a) const
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
template <typename T>
Matrix<T> Matrix<T>::operator/(const int a) const
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
template <typename T>
bool Matrix<T>::operator==(const Matrix& other) const 
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
template <typename T>
bool Matrix<T>::operator==(int scalar) const 
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
template <typename T>
void Matrix<T>::operator=(const Matrix& m)
{
    this.N = m.get_rows();
    this.M = m.get_cols();
    this.data = m.get_data();
}

// Returns inverted matrix.
template <typename T>
Matrix<T> Matrix<T>::operator!()
{
    double det = get_det();
    if (!det) throw std::runtime_error("Matrix cannot be reversed.");

    Matrix newm(N, M);
    transpose();

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
            newm.set(i, j, get_minor(i, j));
    }

    transpose();
    return newm / det;
}

// Print matrix
template <typename T>
std::ostream& operator<< (std::ostream &out, const Matrix<T> &m)
{
    for (int i = 0; i < m.get_rows(); ++i)
    {
        for (int j = 0; j < m.get_cols(); ++j)
            out << m.get_item(i, j) << " ";
        out << std::endl;
    }
    return out;
}

// Create identity matrix
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

int main() {
    Matrix<int> a("data.txt");
    
    std::cout << a.get_det() << std::endl;

    return 0;
}