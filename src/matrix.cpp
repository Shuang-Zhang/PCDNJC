/**
 * \file      matrix.cpp
 * \author    whu.math
 * \copyright WHU resvered
 * \version   1.1
 * \date      2020-3-7
 *
 * \brief     这个矩阵类的内部数据将以列优先，用一维数组存储。目前该类只能实例化float及double类型，后续可以再进行拓展
 * \details
 *
 * \see       其他相关类、函数、变量，文献、书籍、网页等参考资料；
 * \note      依赖关系(用到的第三方库，项目里的其他函数/类)；
 *            测试脚本需要说明测试条件：
 *            操作系统，开发环境，硬件环境，测试数据集，测试结果；
 *            其他不适合放在details里的事项；
 * \bug       已知bug；
 * \todo      未完成事项；
 */
#include "matrix.h"
#include "openblas/cblas.h"
#include "openblas/lapacke.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string.h>
#include <time.h>

#include"utils.h"
    Check::Check() { this->flag = 0; }

    int Check::CheckFlag() const { return this->flag; }

    const Check &Check::SetFlag(int flag, std::string str) const {
        this->strError = str;
        std::cout << this->strError << std::endl;
        this->flag = flag;

        return *this;
    }

    Matrix::Matrix() {
        this->shape = {0, 0};
        this->size  = 0;
        this->pData = nullptr;
    }

    
    Matrix::~Matrix() {
        FreeSpace();
    }

    
    bool Matrix::AllocSpace() {
        if (size == 0) {
            return 0;
        }
        int cols = this->shape[1];
        // pData=new double[size];
        this->pData = new (std::nothrow) double[size];
        double * temp = pData;
        this->vecpData.resize(this->shape[0]);
        for(int i=0;i<this->shape[0];++i){
            this->vecpData[i] = temp;
            temp += cols; 
        }
        if (pData == nullptr) {
            return 1;
        }
        return 0;
    }

    
    bool Matrix::FreeSpace() {
        if (pData != nullptr) {
            delete[] pData;
            pData = nullptr;
        }
        this->vecpData.clear();
        shape = {0, 0};
        size  = 0;
        return 0;
    }


    Matrix::Matrix(const Matrix &inputArray) {
        Array(inputArray);
    }
    Matrix::Matrix(const std::vector<std::vector<double>> &inputArray){
        Array(inputArray);
    }
    Matrix &Matrix::Array(double *inputArray, std::vector<int> shape){
        FreeSpace();
        this->shape = shape;
        size        = 1;
        for (unsigned int i = 0; i < shape.size(); ++i) {
            size *= shape[i];
        }
        if (!AllocSpace()) {
            cblas_dcopy(size, inputArray, 1, pData, 1);
        }
        return *this;




    }
    Matrix &Matrix::Array(const std::vector<std::vector<double>> &inputArray) {
        FreeSpace();
        if (inputArray.empty() || inputArray[0].empty()) {
            this->size = 0;
            return *this;
        }
        this->size = 1;
        this->size *= inputArray[0].size();
        this->size *= inputArray.size();
        this->shape = {(int)inputArray.size(), (int)inputArray[0].size()};
        if (!AllocSpace()) {
            for (int i = 0; i < this->shape[0]; ++i) {
                cblas_dcopy(this->shape[1], inputArray[i].data(), 1, pData + i * this->shape[1], 1);
            }
            
        }
    
        return *this;
    }

    Matrix &Matrix::Array(const Matrix &inputArray) {
        if (&inputArray == this) {
            return *this;
        }
        Empty(inputArray.shape);
        if (size == 0) {
            return *this;
        }
        cblas_dcopy(size, inputArray.pData, 1, pData, 1);
        return *this;
    }

    Matrix &Matrix::operator=(const Matrix &inputArray) {
        return Array(inputArray);
    }
    
    Matrix &Matrix::Zeros(const std::vector<int>& shape){
        Empty(shape);
        for(int i=0;i<this->size;++i){
            this->pData[i] = 0;
        }
        return *this;
    }
    Matrix &Matrix::Empty(const std::vector<int>& shape) {
        FreeSpace();
        this->shape = shape;
        size        = 1;
        for (unsigned int i = 0; i < shape.size(); ++i) {
            size *= shape[i];
        }
        if (!AllocSpace()) {
            return *this;
        }
        return *this;
    }

    int Matrix::Size(int axis) const {
        if (axis == -1) {
            return size;
        }
        if (axis == 0) {
            return this->shape[0];
        }
        if (axis == 1) {
            return this->shape[1];
        }
        return size;
    }

    Matrix &Matrix::Add(const Matrix &matrix) {
        if (shape[0] == matrix.shape[0] && shape[1] == matrix.shape[1]) {
            cblas_daxpy(size,1., matrix.pData, 1, pData, 1);
        }           //要加的矩阵只有一列，原矩阵各列与其相减
        else if (shape[0] == matrix.shape[0] && matrix.shape[1] == 1) {
                for (int i = 0; i < shape[0]; ++i) {
                    for (int j = 0; j < shape[1]; ++j) {
                        pData[i * shape[1] + j] += matrix.pData[i];
                    }
                }
        }
        //要加的矩阵只有一行，原矩阵各行与其相减
        else if (shape[1] == matrix.shape[1] && matrix.shape[0] == 1) {
            for (int i = 0; i < shape[0]; ++i) {
                for (int j = 0; j < shape[1]; ++j) {
                    pData[i * shape[1] + j] += matrix.pData[j];
                }
            }
        } else {
            this->SetFlag(2, "错误:矩阵相加大小不合");
        }
        return *this;
    }

    
    Matrix &Matrix::Add(double b) {
        for (int i = 0; i < size; ++i) {
            pData[i] += b;
        }
        return *this;
    }

    
    Matrix &Matrix::Subtract(const Matrix &matrix) {
        if (shape[0] == matrix.shape[0] && shape[1] == matrix.shape[1]) {
            cblas_daxpy(size,-1., matrix.pData, 1, pData, 1);
        }           //要加的矩阵只有一列，原矩阵各列与其相减
        else if (shape[0] == matrix.shape[0] && matrix.shape[1] == 1) {
                for (int i = 0; i < shape[0]; ++i) {
                    for (int j = 0; j < shape[1]; ++j) {
                        pData[i * shape[1] + j] -= matrix.pData[i];
                    }
                }
        }
        //要加的矩阵只有一行，原矩阵各行与其相减
        else if (shape[1] == matrix.shape[1] && matrix.shape[0] == 1) {
            for (int i = 0; i < shape[0]; ++i) {
                for (int j = 0; j < shape[1]; ++j) {
                    pData[i * shape[1] + j] -= matrix.pData[j];
                }
            }
        } else {
            this->SetFlag(2, "错误:矩阵相加大小不合");
        }
        return *this;
    }

    
    Matrix &Matrix::Subtract(double b) {
        for (int i = 0; i < size; ++i) {
            pData[i] -= b;
        }
        return *this;
    }

    
    Matrix &Matrix::Multiply(double b) {
        cblas_dscal(size, b, pData, 1);
        return *this;
    }

    
    Matrix &Matrix::Divide(double b) {
        double temp = 1. / b;
        cblas_dscal(size, temp, pData, 1);
        return *this;
    }

    Matrix &Matrix::Dot(const Matrix &matrix) {
        if (shape[0] == matrix.shape[0] && shape[1] == matrix.shape[1]) {
            for (int i = 0; i < size; ++i) {
                pData[i] *= matrix.pData[i];
            }
        }else if (shape[0] == matrix.shape[0] && matrix.shape[1] == 1) {
            for (int i = 0; i < shape[0]; ++i) {
                for (int j = 0; j < shape[1]; ++j) {
                    pData[i * shape[1] + j] *= matrix.pData[i];
                }
            }
        }
        //原矩阵各行分别与新矩阵做点乘
        else if (shape[1] == matrix.shape[1] && matrix.shape[0] == 1) {
            for (int i = 0; i < shape[0]; ++i) {
                for (int j = 0; j < shape[1]; ++j) {
                    pData[i * shape[1] + j] *= matrix.pData[j];
                }
            }
        } else {
            this->SetFlag(2,"错误:矩阵点乘大小不合");
            return *this;
        }
        return *this;
    }

    
    Matrix &Matrix::Multiply(const Matrix &matrix) {
        if (this->shape[1] != matrix.shape[0]) {
            std::cout<<"错误: 矩阵相乘大小不合"<<std::endl;;
            return *this;
        }
        Matrix temp(*this);
        Zeros({this->shape[0], matrix.shape[1]});{
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, temp.shape[0], matrix.shape[1], temp.shape[1], 1,
                       temp.pData, temp.shape[1], matrix.pData, matrix.shape[1], 0, pData, this->shape[1]);
        } 
        return *this;
    }

    
    Matrix &Matrix::Divide(const Matrix &matrix) {
        if (matrix.shape[0] < matrix.shape[1]) {
            this->SetFlag(2, "错误:系数方程组行数过低");
            return *this;
        }
        //运行LAPACK_gels 会改变被除矩阵 matrix数据，所以需要对matrix进行拷贝
        //另外，虽然其将计算结果是直接保存在原矩阵上，
        //但原矩阵大小可能与被除结果不一致，所以还是要重新分配大小，进行拷贝
        Matrix temp(*this);
        Matrix tempMatrix(matrix);
        Zeros({tempMatrix.shape[1], this->shape[1]});
        LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', tempMatrix.shape[0], tempMatrix.shape[1], temp.shape[1],
                        tempMatrix.pData, tempMatrix.shape[1], temp.pData, temp.shape[1]);
        cblas_dcopy(this->size, temp.pData, 1, pData, 1);
        return *this;
    }

    Matrix &Matrix::Transpose() {
        if (this->size != 0) {
            Matrix temp(*this);
            if (this->shape[0] == this->shape[1]) {
                cblas_domatcopy(CblasRowMajor, CblasTrans, this->shape[1], this->shape[0], 1, temp.pData,
                                this->shape[1], pData, this->shape[1]);
            }
            else {
                if (this->shape[0] != 1 && this->shape[1] != 1) {
                    int startIndex = 0;
                    for (int i = 0; i < this->shape[1]; ++i) {
                        cblas_dcopy(this->shape[0], temp.pData + i, this->shape[1], pData + startIndex, 1);
                        startIndex += shape[0];
                    }
                }
                this->shape[0] = temp.shape[1];
                this->shape[1] = temp.shape[0];
            }
        }
        return *this;
    }

    Matrix &Matrix::Get(Matrix &result, const std::vector<int> &cols, int axis) const {
        int ilength = (int)cols.size();
        if (axis == 0) {
            result.Empty({shape[0], ilength});
            for (int i = 0; i < ilength; ++i) {
                if (cols[i] < 0 || cols[i] >= Size(1)) {
                    result.SetFlag(2, "错误：矩阵按列截取超出索引");
                    return result;
                }
                cblas_dcopy(shape[0], pData + cols[i], shape[1], result.pData + i, ilength);
                
            }
        } else if (axis == 1) {
            result.Empty({ilength, shape[1]});
            for (int i = 0; i < ilength; ++i) {
                if (cols[i] < 0 || cols[i] >= Size(0)) {
                    result.SetFlag(2, "错误：矩阵按行截取超出索引");
                    return result;
                }
                cblas_dcopy(shape[1], pData + cols[i] * shape[1], 1, result.pData + i * shape[1], 1);
                
            }
        } else {
            result.SetFlag(2, "错误：矩阵截取axis 只能取 0,1");
            return result;
        }
        return result;
    }
    Matrix &Matrix::Get(Matrix &result, int begin, int end, int axis) const {
        if (begin > end) {
            result.SetFlag(2, "begin > end,截取矩阵错误");
            return result;
        }
        std::vector<int> index;
        for (int i = begin; i < end; ++i) {
            index.push_back(i);
        }
        Get(result, index, axis);
        return result;
    }

    
    std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
        int irows = matrix.Size(0);
        int icols = matrix.Size(1);
        for (int i = 0; i < irows; ++i) {
            int j=0;
            for (; j < icols-1; ++j) {
                os << matrix[i][j] << " ";
            }
            if(j<icols){
                os <<matrix[i][j];
            }
            os << std::endl;
        }
        

        return os;
    }

    std::ofstream &operator<<(std::ofstream & os, const Matrix &matrix){
        int irows = matrix.Size(0);
        int icols = matrix.Size(1);
        for (int i = 0; i < irows; ++i) {
            for (int j = 0; j < icols; ++j) {
                os << matrix[i][j] << " ";
            }
            os << std::endl;
        }
        return os;
    }
bool IsEqual(const Matrix &matrix1, const Matrix &matrix2, double epsilon){
        if ((matrix1.Size(0) != matrix2.Size(0)) || (matrix1.Size(1) != matrix2.Size(1))) {
            std::cout << "输入的两个矩阵不相等,分别为" << std::endl;
            std::cout << "matrix1:" << std::endl << matrix1 << "matrix2:" << std::endl << matrix2 << std::endl;
            return false;
        }
        if (matrix1.Size() == 0) {
            return true;
        }
        double *p1 = matrix1.GetPtr();
        double *p2 = matrix2.GetPtr();
        for (int i = 0; i < matrix1.Size(); ++i) {
            if (std::fabs(*(p1 + i) - *(p2 + i)) < epsilon) {
            } else {
                if ((IsNAN(*(p1 + i)) && IsNAN(*(p2 + i))) || (IsInf(*(p1 + i)) && IsInf(*(p2 + i)))) {
                    continue;
                }
                std::cout << "输入的两个矩阵不相等，分别为" << std::endl;
                std::cout << "matrix1:" << std::endl << matrix1 << "matrix2:" << std::endl << matrix2 << std::endl;
                return false;
            }
        }
        return true;
}
    
    Matrix &Matrix::Sum(Matrix &sum, int axis) const {
        if (axis == -1) {
            sum.Zeros({1, 1});
            for (int j = 0; j < size; ++j) {
                sum.pData[0] += pData[j];
            }
        } else if (axis == 0) {
            sum.Zeros({1, this->shape[1]});
            for (int i = 0; i < this->shape[0]; ++i) {
                for (int j = 0; j < this->shape[1]; ++j) {
                    sum.pData[j] += pData[i * this->shape[1] + j];
                }
            }
            
        } else if (axis == 1) {
            sum.Zeros({this->shape[0], 1});
            for (int i = 0; i < this->shape[0]; ++i) {
                for (int j = 0; j < this->shape[1]; ++j) {
                    sum.pData[i] += pData[i * shape[1] + j];
                }
            }
            
        } else {
            sum.SetFlag(2, "sum 求和 axis 只能设置为 -1，0, 1");
        }
        return sum;
    }
    
    Matrix &Matrix::Max(Matrix &max, int axis) const {
        std::vector<int> argmax;
        // true 表示求最大值
        Argmax(max, argmax, axis, true);
        return max;
    }
    /**
     * @brief 求最大值元素的下标
     */
    
    std::vector<int> &Matrix::Argmax(std::vector<int> &argmax, int axis) const {
        Matrix max;
        argmax.clear();
        // true 表示求最大值
        Argmax(max, argmax, axis, true);
        return argmax;
    }

    
    Matrix &Matrix::Mean(Matrix &mean, int axis) const {
        Sum(mean, axis);
        if (axis == -1) {
            mean.Divide((double)size);
        } else if (axis == 0) {
            mean.Divide((double)shape[0]);
        } else if (axis == 1) {
            mean.Divide((double)shape[1]);
        } else {
            mean.SetFlag(2, "求mean,axis只能输入-1 或 0 或1");
        }

        return mean;
    }

    
    Matrix &Matrix::Average(Matrix &average, int axis, const Matrix &sampleWeight) const {
        Matrix sum;
        Matrix x(*this);
        average.FreeSpace();
        // axis = -1 且 sampleWeight 非空的情况尚未考虑
        if (sampleWeight.Size() != 0) {
            sampleWeight.Sum(sum, -1);
            Matrix sampleWeightCopy(sampleWeight);
            sampleWeightCopy.Divide(sum.pData[0]);
            if (axis == 0) {
                sampleWeightCopy.Reshape({x.Size(0), 1});
            } else if (axis == 1) {
                sampleWeightCopy.Reshape({1, x.Size(1)});
            }

            x.Dot(sampleWeightCopy);
            x.Sum(average, axis);
        } else {
            x.Mean(average, axis);
        }
        return average;
    }

    
    Matrix &Matrix::Var(Matrix &var, int axis, const Matrix &mean, int ddof) const {
        Matrix mu;
        if (mean.Size() == 0) {
            Mean(mu, axis);
        } else {
            mu.Array(mean);
        }
        // axis =-1情况不用区分行优先列优先
        if (axis == -1) {
            var.Zeros({1, 1});
            double fVar = 0;
            for (int i = 0; i < size; ++i) {
                fVar += pow((pData[i] - mu.pData[0]), 2);
            }
            fVar /= size - ddof;
            var.pData[0] = fVar;
        } else if (axis == 0) {
            var.Zeros({1, shape[1]});
            for (int i = 0; i < shape[0]; ++i) {
                int startIndex = i * shape[1];
                for (int j = 0; j < shape[1]; ++j) {
                    var.pData[j] += pow(pData[j + startIndex] - mu.pData[j], 2);
                }
            }
            var.Divide(shape[0] - ddof);
        } else if (axis == 1) {
            var.Zeros({shape[0], 1});
            for (int i = 0; i < shape[0]; ++i) {
                int startIndex = i * shape[1];
                //保证内层循环依然是连续读写数据
                for (int j = 0; j < shape[1]; ++j) {
                    var.pData[i] += pow(pData[startIndex + j] - mu.pData[i], 2);
                }
            }
            var.Divide(shape[1] - ddof);
        } else {
            var.SetFlag(2, "axis 输入值只能为 -1,0,1");
        }
        return var;
    }

    
    Matrix &Matrix::Std(Matrix &var, int axis, const Matrix &mean) const {
        Var(var, axis, mean);
        for (int i = 0; i < var.Size(); ++i) {
            *(var.pData + i) = sqrt((*(var.pData + i)));
        }
        return var;
    }

    Matrix &Matrix::Linspace(double start, double stop, int num, bool endpoint) {
        Empty({1, num});
        float step = 0;
        if (endpoint == 1) {
            step = (stop - start) / (num - 1);
        } else {
            step = (stop - start) / num;
        }
        for (int i = 0; i < size; ++i) {
            pData[i] = start + i * step;
        }
        return *this;
    }
    Matrix &Matrix::Exp(){
        for(int i=0;i<this->size;++i){
            this->pData[i] = exp(this->pData[i]);
        }
        return *this;
    }
    //后续改进成按行去重或按列去重
    Matrix &Matrix::MyUnique(Matrix &result) const {
        //仅针对一维数组使用该函数
        Matrix temp(*this);
        std::sort(temp.pData, temp.pData + size);
        // unique(a,a+n)返回的是从去重后（不重复数列中）最后一个元素的尾地址，减去地址a就是数列的长度
        int length = std::unique(temp.pData, temp.pData + size) - temp.pData;
        std::vector<int> shape;
        if (this->shape[0] == 1) {
            shape = {1, length};
        } else {
            shape = {length, 1};
        }

        result.Array(temp.pData, shape);

        return result;
    }
     
    Matrix &Matrix::Concatenate(Matrix &result, const Matrix &matrix, int axis) const {
        if (this->Size() != 0 && matrix.Size() != 0) {
            if (axis == 0) {
                if (Size(1) != matrix.Size(1)) {
                    result.SetFlag(2, "失败:矩阵拼接列数不一致");
                    return result;
                }
                result.Empty({this->shape[0] + matrix.shape[0], this->shape[1]});
                    //行优先矩阵对新矩阵拼接在原矩阵下方速度比较快
                    cblas_dcopy(this->size, pData, 1, result.pData, 1);
                    cblas_dcopy(matrix.size, matrix.pData, 1, result.pData + this->size, 1);

            } else if (axis == 1) {
                if (Size(0) != matrix.Size(0)) {
                    result.SetFlag(2, "失败:矩阵拼接行数不一致");
                    return result;
                }
                result.Empty({this->shape[0], this->shape[1] + matrix.shape[1]});

                //原矩阵每次拷贝起始位置
                int indexThis = 0;
                // matrix拷贝起始位置
                int indexMatrix = 0;
                // result 拷贝起始位置
                int indexResult = 0;
                //一行一行进行复制
                for (int i = 0; i < this->shape[0]; ++i) {
                    cblas_dcopy(this->shape[1], pData + indexThis, 1, result.pData + indexResult, 1);
                    indexResult += this->shape[1];
                    cblas_dcopy(matrix.shape[1], matrix.pData + indexMatrix, 1, result.pData + indexResult, 1);
                    indexThis += this->shape[1];
                    indexMatrix += matrix.shape[1];
                    indexResult += matrix.shape[1];
                }
                
            }
        } else if (this->size == 0) {
            result.Array(matrix);
        }
        return result;
    }

    Matrix &Matrix::Concatenate(const Matrix &matrix,int axis =0){
        Matrix temp(*this);
        return temp.Concatenate(*this,matrix,axis);       
    }

    
    Matrix &Matrix::Reshape(std::vector<int> shape) {
        int flag      = -1;
        int zerocount = 0;
        int newsize   = 1;
        for (unsigned int i = 0; i < shape.size(); ++i) {
            if (shape[i] != -1) {
                newsize *= shape[i];
            } else {
                flag = i;
                ++zerocount;
            }
        }
        if (zerocount > 1) {
            this->SetFlag(2, "reshape error,输入的shape中-1的个数超过1");
            return *this;
        }
        if (zerocount == 0 && newsize != this->size) {
            this->SetFlag(2, "reshape size error:新的shape的大小与之前的不一致");
            return *this;
        }
        if (zerocount == 1) {
            shape[flag] = this->size / newsize;
            this->shape = shape;
        } else if (zerocount == 0) {
            this->shape = shape;
        }
        return *this;
    }
    Matrix &Matrix::Min(Matrix &min, int axis) const {
        std::vector<int> argmin;
        Argmax(min, argmin, axis, false);
        return min;
    }

    // n次抽取随机数，每次从[0,n-i-1) 选择一个数 与 第 n-i-1个数交换，i=0,1,...n-1
    std::vector<int> &Matrix::RandVector(std::vector<int> &index, int n) {
        index.clear();
        index.reserve(n);
        for (int i = 0; i < n; ++i) {
            index.push_back(i);
        }
        for (int i = n; i > 0; --i) {
            int id = (int)(rand() % i);
            //交换值
            int value    = index[i - 1];
            index[i - 1] = index[id];
            index[id]    = value;
        }
        return index;
    }
    
    Matrix &Matrix::Rand(std::vector<int> shape) {
        this->Empty(shape);

        double *p = GetPtr();
        for (int i = 0; i < this->size; ++i) {
            *(p + i) = (double)rand();
        }
        return *this;
    }
    bool Matrix::SaveTXT(std::string strfile, std::string strmatrix) {
        std::ofstream outfile;
        outfile.open(strfile);
        outfile << strmatrix << "：" << this->shape[0] << "," << this->shape[1] << std::endl;
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                outfile << (*this)[i][j] << " ";
            }
            outfile << std::endl;
        }
        outfile.close();
        return 0;
    }
    std::vector<int> & Matrix::Matrix2VectorInt(std::vector<int> &result){
        result.clear();
        for (int i = 0; i < this->size;++i) {
            result.push_back(pData[i]);
        }
        return result;        
    }
    
    std::vector<std::vector<double>> &Matrix::Matrix2vector2(std::vector<std::vector<double>> &result) const {
        result.clear();
        std::vector<double> temp;
        int startIndex = 0;
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                temp.emplace_back(*(pData + startIndex + j));
            }
            startIndex += shape[1];
            result.emplace_back(std::move(temp));
        }
        return result;
    }

    Matrix & Matrix::RandN(const std::vector<int>& shape,int seednum){
        Empty(shape);
        srand(seednum);
        for(int i=0;i<this->size;++i){static double V1, V2, S;
            static int phase = 0;
            double X;
            
            if ( phase == 0 ) {
                do {
                    double U1 = (double)rand() / RAND_MAX;
                    double U2 = (double)rand() / RAND_MAX;
                    
                    V1 = 2 * U1 - 1;
                    V2 = 2 * U2 - 1;
                    S = V1 * V1 + V2 * V2;
                } while(S >= 1 || S == 0);
                
                X = V1 * sqrt(-2 * log(S) / S);
            } else
                X = V2 * sqrt(-2 * log(S) / S);
                
            phase = 1 - phase;
            pData[i] = X;
        }
        return *this;
    }
 
    Matrix & Matrix::Chol(Matrix &result) const{
        result.Array(*this);
        char uplo = 'L';
        int order = result.Size(0);
        int lda = order;
        int info = 0;
        LAPACK_dpotrf(&uplo,&order,result.GetPtr(),&lda,&info);
            //将下三角元素置为 0
        double * row = result.GetPtr();
        for(int i=0;i<order;++i){
            for(int j=0;j<i;++j){
                row[j] = 0;
            }
            row += order;
        }
        return result;


    }
    double *Matrix::GetPtr() const{
        return this->pData;
    }

    Matrix &Matrix::Argmax(Matrix &max, std::vector<int> &argmax, int axis, bool status) const {
        if (this->size == 0) {
            max.Empty({0, 0});
            return max;
        }
        //如果是对整个矩阵求最大值元素及其下标
        if (axis == -1) {
            max.Empty({1, 1});
            double tmax      = pData[0];
            int targmax = 0;
            for (int i = 0; i < size; ++i) {
                if (MyCompare(pData[i], tmax, status)) {
                    tmax    = pData[i];
                    targmax = i;
                }
            }
            max.pData[0] = tmax;
            argmax.push_back(targmax);
        }else  if (axis == 0) {
            max.Empty({1, shape[1]});
            for (int i = 0; i < shape[1]; ++i) {
                max.pData[i] = pData[i];
            }
            argmax.resize(shape[1]);
            for (int i = 1; i < shape[0]; ++i) {
                int startIndex = i * shape[1];
                //内循环保证连续读取
                for (int j = 0; j < shape[1]; ++j) {
                    if (MyCompare(pData[startIndex + j], max.pData[j], status)) {
                        max.pData[j] = pData[startIndex + j];
                        argmax[j]    = i;
                    }
                }
            }
        }
        //如果是对矩阵各行求最大值
        else if (axis == 1) {
            max.Empty({shape[0], 1});
            for (int i = 0; i < shape[0]; ++i) {
                int startIndex = i * shape[1];
                int tmax       = pData[startIndex];
                int targmax    = 0;
                //内循环保证元素连续读取
                for (int j = 1; j < shape[1]; ++j) {
                    if (MyCompare(pData[startIndex + j], tmax, status)) {
                        tmax    = pData[startIndex + j];
                        targmax = j;
                    }
                }
                max.pData[i] = tmax;
                argmax.push_back(targmax);
            }
        } else {
            max.SetFlag(2, "axis 输入值为 -1,0,1");
        }
        return max;
    }
    
 
