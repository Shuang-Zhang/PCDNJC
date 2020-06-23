
#ifndef MATRIX
#define MATRIX

#include <fstream>
#include <ostream>
#include <iostream>

#include <string>
#include <vector>
#include<math.h>
    /**
     * \brief 处理错误信息，适时进行返回
     */
class Check {
    public:
    Check();
    /**
     * \brief 检查当前类的状态，一般flag为0表示正常，1表示警告,2表示错误,也可以根据自身需要再确定flag的实际含义
     * \return 返回 this->flag
     */
    int CheckFlag() const;

    /**
     * \brief  修改类的状态
     * \param flag flag 的值为0时表示正常，1表示警告,2表示错误
     * \param str 可将错误信息存入，若其非空，则会在屏幕上打印出该错误信息
     * \return 返回 *this
     */
    const Check &SetFlag(int flag, std::string str = "") const;

    protected:
    // 判断当前矩阵状态是否正常运行出来的结果的变量
    // flag 的值为0时表示正常，1表示警告,2表示错误
    mutable int flag;
    // 可将相关错误信息储存在该变量
    mutable std::string strError;
};
// PI
const double ML_PI = 3.1415926535;

class Matrix:public Check{
    public:
    Matrix();
    ~Matrix();
    Matrix(const std::vector<std::vector<double>> &inputArray);
    Matrix(const std::vector<double> &inputArray);
    Matrix(const Matrix &inputArray);
    Matrix &Array(double *inputArray, std::vector<int> shape);
    Matrix &Array(const std::vector<std::vector<double>> &inputArray);
    Matrix &Array(const Matrix &inputArray);
    Matrix &operator=(const Matrix &);
    Matrix &Zeros(const std::vector<int>& shape);
    Matrix &Empty(const std::vector<int>& shape);
    int Size(int axis = -1) const;

    // 自增运算
    Matrix &Add(const Matrix &matrix);
    Matrix &Add(double b);

    // 自减运算
    Matrix &Subtract(const Matrix &matrix);
    Matrix &Subtract(double b);

    // 自乘运算
    Matrix &Multiply(double b);

    // 自除运算
    Matrix &Divide(double b);

    // 矩阵点乘
    Matrix &Dot(const Matrix &matrix);
    // 矩阵相乘,相乘的结果存储在 result中
    Matrix &Multiply(const Matrix &matrix);
    // 矩阵相除，记原矩阵为b，即  A*x=b ，求 x=b/A ,原矩阵变为x
    Matrix &Divide(const Matrix &A);
    /**
     * \brief 将原矩阵转置
     */
    Matrix &Transpose();
    /**
     * \brief 矩阵元素求和
     * \param axis用于二维矩阵时，axis = -1表示整个矩阵的求和，axis = 0表示原矩阵表示按列求和，axis =
     * 1原矩阵表示按行求和
     */
    Matrix &Sum(Matrix &sum, int axis = -1) const;

    /**
     * \brief 矩阵元素求最大值
     * \param axis 用于二维矩阵时，axis = -1 求整个矩阵元素的最大值，axis = 0表示原矩阵按列求最大值，axis =
     * 1原矩阵按行求最大值
     */
    Matrix &Max(Matrix &max, int axis = -1) const;

    /**
     * \brief 求最大值元素的下标
     */
    std::vector<int> &Argmax(std::vector<int> &argmax, int axis = -1) const;
    
    /**
     * \brief 求矩阵元素的均值
     * \param axis 用于二维矩阵时， axis = -1 ，axis =  0 表示原矩阵按列求均值，axis = 1 表示原矩阵按行求均值
     */
    Matrix &Mean(Matrix &mean, int axis = -1) const;

    /**
     * \brief 求矩阵元素的均值，与Mean 类似，主要为样本加权，首先会对每个样本即每行乘上其对应权重，然后再求平均
     */
    Matrix &Average(Matrix &average, int axis = -1, const Matrix &sampleWeight = Matrix()) const;

    /**
     * \brief 求矩阵元素的方差
     * \param axis 用于二维矩阵时 axis = -1 表示求整个矩阵元素的方差，axis =  0表示原矩阵按列求方差，axis =
     *             1表示原矩阵按行求方差 ddof 表示delta自由度
     * \param mean 可以将平均值作为参数传入，将减少再次计算均值
     * \param ddof 是“Delta Degrees of
     * Freedom”，表示自由度的个数，在计算方差时，分子是各个值和均值的差的平方之和，分母为（N-ddof）
     */
    Matrix &Var(Matrix &result, int axis = -1, const Matrix &mean = Matrix(), int ddof = 0) const;

    // sqrt(var)
    Matrix &Std(Matrix &result, int axis = -1, const Matrix &mean = Matrix()) const;

    Matrix &Linspace(double start, double stop, int num = 50, bool endpoint = 1);
    
    Matrix &Exp();
    /**
     * \brief 返回去重后的matrix先仅用于一维向量，会去掉矩阵内部重复的元素
     */
    Matrix &MyUnique(Matrix &result) const;

    /**
     * \brief 数组拼接，axis = 0 是矩阵拼接在原矩阵下方，axis = 1 是矩阵拼接在原矩阵右边
     */
    Matrix &Concatenate(Matrix &result, const Matrix &matrix, int axis = 0) const;
    Matrix &Concatenate(const Matrix &matrix,int axis);
    /**
     * \brief 重新设置矩阵形状
     */
    Matrix &Reshape(std::vector<int> shape);
    inline double &Get(const int& i) const{
       return pData[i];
    }
    Matrix &Get(Matrix &result, int begin, int end, int axis) const;
    
    Matrix &Get(Matrix &result, const std::vector<int> &cols, int axis = 0) const;
    /**
     * \param axis 当 axis 为-1，求取所有元素最小值
             当axis为0,求取各列最小值,返回大小为 (1,this->shape[1])，
                当axis为1,求取各行最小值,返回大小为 (this->shape[0],1)
        */
    Matrix &Min(Matrix &result, int axis = -1) const;

    /**
     * \brief 生成指定大小的均匀分布的的矩阵
     */
    Matrix &Rand(std::vector<int> shape);
    /**
     * \brief 生成一组从0到n-1的随机序列，长度为n，且各数据不重复，可用于将一个矩阵按行或者按列打乱
     */
    static std::vector<int> &RandVector(std::vector<int> &newIndex, int n);

    /**
     * \brief 将矩阵存为文本，但不会存下原来矩阵order
     * \param strmatrix 可以存下被存矩阵的名字，以后可以按名字从从文本读出对应矩阵
     */
    bool SaveTXT(std::string strfile, std::string strmatrix = "matrix");
    /**
     *  \brief 将原矩阵转换为 vector<vector<double>> 容器存取
     */
    std::vector<std::vector<double>> &Matrix2vector2(std::vector<std::vector<double>> &result) const;
    
    std::vector<int> & Matrix2VectorInt(std::vector<int> &result);
    
    Matrix & RandN(const std::vector<int>& shape,int seednum =0); 

    Matrix & Chol(Matrix &result) const;
    double *GetPtr() const;
    inline double* operator[](const int& i) const{
        return this->vecpData[i];
    }
    inline double* operator[](const int& i){
        return this->vecpData[i];
    }
    private:
    // 申请分配内存空间
    bool AllocSpace();

    // 释放内存空间
    bool FreeSpace();

    Matrix &Argmax(Matrix &max, std::vector<int> &argmax, int axis, bool status) const;

    // 用于结合Max 和 Min 代码
    inline bool MyCompare(const double &a, const double &b, bool status) const {
        if (status == true) {
            return a > b;
        } else {
            return a < b;
        }
    }
    // 数组头指针
    double *pData = nullptr;
    std::vector<double *> vecpData;
    // 用于表示矩阵大小,shape[i]是第i维度的大小
    // 用于二维矩阵时，shape[0]是行数，m_shape[1]是列数
    std::vector<int> shape;

    // 矩阵的大小
    int size;
};

// 作用:Matrix 内元素的输入和打印
// 使用方法 : cout<< Matrix对象;    //可打印出Matrix里的元素
std::ostream &operator<<(std::ostream &, const Matrix &);
std::ofstream &operator<<(std::ofstream &, const Matrix &);
template <typename T>
bool IsNAN(T x) {
    return x != x;
}

template <typename T>
bool IsInf(T x){
        return !IsNAN(x) && IsNAN(x - x);
}
/**
 * \brief 判断两个矩阵数据是否相同 ,矩阵对应元素值相差epsilon以内认为是相等的
 */
bool IsEqual(const Matrix &matrix1, const Matrix &matrix2, double epsilon = 1e-3);
template<typename T>
bool IsEqual(const std::vector<T> vec1, const std::vector<T> vec2, double epsilon = 1e-3){
    if (vec1.size() != vec2.size()) {
        std::cout << "输入的两个vector不相等，分别为" << std::endl;
        std::cout << "vec1:" << std::endl << vec1 << "vec2:" << std::endl << vec2 << std::endl;
        return false;
    }
    for (unsigned int i = 0; i < vec1.size(); ++i) {
        if (fabs(vec1[i] - vec2[i]) > epsilon) {
            std::cout << "输入的两个vector不相等，分别为" << std::endl;
            std::cout << "vec1:" << std::endl << vec1 << "vec2:" << std::endl << vec2 << std::endl;
            return false;
        }
    }
    return true;
}

// 打印vector值，用于调试
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &m) {
    for (unsigned int i = 0; i < m.size(); ++i) {
        os << m[i] << " ";
    }
    std::cout << std::endl;
    return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<std::vector<T>> &m) {
    for (unsigned int i = 0; i < m.size(); ++i) {
        os << m[i];
    }
    std::cout << std::endl;
    return os;
}


#endif
