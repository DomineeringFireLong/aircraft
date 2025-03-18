#include "Basic.h"

#define IS_ELLIPSOID true     // 是否考虑地球扁率
/***********************************************************************
	BDPoint3
	**********************************************************************/

BDPoint3::BDPoint3()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

BDPoint3::BDPoint3(const BDPoint3& point)
{
	x = point.x;
	y = point.y;
	z = point.z;
}

BDPoint3::BDPoint3(double dX, double dY, double dZ)
{
	x = dX;
	y = dY;
	z = dZ;
}

BDPoint3::BDPoint3(const double data[3])
{
	x = data[0];
	y = data[1];
	z = data[2];
}

BDPoint3::BDPoint3(const std::vector<double>& i)
{
	assert(i.size() >= 3);
	x = i[0];
	y = i[1];
	z = i[2];
}

BDPoint3::~BDPoint3()
{}

BDPoint3& BDPoint3::operator=(const BDPoint3& point)
{
	x = point.x;
	y = point.y;
	z = point.z;
	return *this;
}

double BDPoint3::Distance(const BDPoint3& point) const
{
	return double(sqrt((point.x - x) * (point.x - x) + (point.y - y) * (point.y - y) + (point.z - z) * (point.z - z)));
}

double	BDPoint3::norm() const
{
	return sqrt(x * x + y * y + z * z);
}

std::string BDPoint3::toStdString() const
{
	return toStdString(x) + " " + toStdString(y) + " " + toStdString(z);
}


BDPoint3& BDPoint3::deg2rad()
{
	x = Deg2Rad(x);
	y = Deg2Rad(y);
	return *this;
}

BDPoint3& BDPoint3::rad2deg()
{
	x = Rad2Deg(x);
	y = Rad2Deg(y);
	return *this;
}

BDPoint3& BDPoint3::alldeg2rad()
{
	x = Deg2Rad(x);
	y = Deg2Rad(y);
	z = Deg2Rad(z);
	return *this;
}

BDPoint3& BDPoint3::allrad2deg()
{
	x = Rad2Deg(x);
	y = Rad2Deg(y);
	z = Rad2Deg(z);
	return *this;
}

void BDPoint3::toArray(double arr[3]) const
{
	arr[0] = x;
	arr[1] = y;
	arr[2] = z;
}

void BDPoint3::fromArray(double arr[3])
{
	this->x = arr[0];
	this->y = arr[1];
	this->z = arr[2];
}

double& BDPoint3::operator[](int iIndex)
{
	switch (iIndex)
	{
	case 0:
		return x;
	case 1:
		return y;
	case 2:
		return z;
	default:
		assert(false);
		return x;
	}
}

double BDPoint3::operator[](int iIndex) const
{
	switch (iIndex)
	{
	case 0:
		return x;
	case 1:
		return y;
	case 2:
		return z;
	default:
		assert(false);
		return 0.0;
	}
}

BDPoint3 BDPoint3::operator+(BDPoint3 point) const
{
	return BDPoint3(x + point.x, y + point.y, z + point.z);
}

BDPoint3 BDPoint3::operator-(BDPoint3 point) const
{
	return BDPoint3(x - point.x, y - point.y, z - point.z);
}

BDPoint3 BDPoint3::operator-() const
{
	BDPoint3 m;
	m.x = -x;
	m.y = -y;
	m.z = -z;
	return m;
}

BDPoint3 BDPoint3::operator*(double dCoe) const
{
	return BDPoint3(x * dCoe, y * dCoe, z * dCoe);
}

double BDPoint3::operator*(BDPoint3 point) const
{
	return double(x * point.x + y * point.y + z * point.z);
}

BDPoint3 BDPoint3::operator/(double a) const
{
	return BDPoint3(x / a, y / a, z / a);
}

bool BDPoint3::operator==(BDPoint3 point) const
{
	if (x == point.x && y == point.y && z == point.z)
		return true;
	else
		return false;
}

bool BDPoint3::operator!=(BDPoint3 point) const
{
	if (x == point.x && y == point.y && z == point.z)
		return false;
	else
		return true;
}

BDPoint3 BDPoint3::cross(const BDPoint3& p) const
{
	BDPoint3 r;
	r[0] = y * p[2] - z * p[1];
	r[1] = z * p[0] - x * p[2];
	r[2] = x * p[1] - y * p[0];
	return r;
}

std::string BDPoint3::toStdString(double d)
{
	std::string str;
	std::stringstream ss;
	ss << d;
	ss >> str;
	return str;
}



/***********************************************************************
BDMatrix
***********************************************************************/

BDMatrix::BDMatrix()
{
	m_nNumColumns = 1;
	m_nNumRows = 1;
	m_pData = nullptr;
	bool bSuccess = Init(m_nNumRows, m_nNumColumns);
	assert(bSuccess);
}

BDMatrix::BDMatrix(unsigned int nRows, unsigned int nCols)
{
	m_nNumRows = nRows;
	m_nNumColumns = nCols;
	m_pData = nullptr;
	bool bSuccess = Init(m_nNumRows, m_nNumColumns);
	assert(bSuccess);
}

BDMatrix::BDMatrix(int nRows, int nCols, double value[])
{
	m_nNumRows = nRows;
	m_nNumColumns = nCols;
	m_pData = nullptr;
	bool bSuccess = Init(m_nNumRows, m_nNumColumns);
	assert(bSuccess);

	SetData(value);
}

// 方阵构造函数
BDMatrix::BDMatrix(int nSize)
{
	m_nNumRows = nSize;
	m_nNumColumns = nSize;
	m_pData = nullptr;
	bool bSuccess = Init(nSize, nSize);
	assert(bSuccess);
}

// 方阵
// 1. int nSize - 方阵行列数;
// 2. double value[] - 一维数组，长度为nRows*nRows，存储方阵各元素的值;
BDMatrix::BDMatrix(int nSize, double value[])
{
	m_nNumRows = nSize;
	m_nNumColumns = nSize;
	m_pData = nullptr;
	bool bSuccess = Init(nSize, nSize);
	assert(bSuccess);

	SetData(value);
}

BDMatrix::BDMatrix(const std::vector<double>& value)
{
	m_pData = nullptr;
	bool bSuccess = Init(1, value.size());
	assert(bSuccess);

	for (size_t i = 0; i < value.size(); i++)
	{
		SetElement(0, i, value[i]);
	}
}

BDMatrix::BDMatrix(const std::vector<std::vector<double> >& value)
{
	m_pData = nullptr;
	if (value.size() == 0) {
		Init(0, 0);
		return;
	}

	bool bSuccess = Init(value.size(), value[0].size());
	assert(bSuccess);

	for (int i = 0; i < GetNumRows(); i++)
	{
		for (int j = 0; j < GetNumColumns(); j++)
			SetElement(i, j, value[i][j]);
	}
}

BDMatrix::BDMatrix(const BDMatrix& other)
{
	m_nNumColumns = other.GetNumColumns();
	m_nNumRows = other.GetNumRows();
	m_pData = nullptr;
	bool bSuccess = Init(m_nNumRows, m_nNumColumns);
	assert(bSuccess);

	// copy the pointer
	memcpy(m_pData, other.m_pData, sizeof(double) * m_nNumColumns * m_nNumRows);
}

BDMatrix::~BDMatrix()
{
	if (m_pData)
	{
		delete[] m_pData;
		m_pData = nullptr;
	}
}

// 初始化函数
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
// 返回值：bool 型，初始化是否成功
bool BDMatrix::Init(int nRows, int nCols)
{
	clear();

	m_nNumRows = nRows;
	m_nNumColumns = nCols;
	int nSize = nCols * nRows;
	if (nSize < 0)
		return false;

	// 分配内存
	m_pData = new double[nSize];

	if (m_pData == nullptr)
		return false;					// 内存分配失败

	// 将各元素值置0
	memset(m_pData, 0, sizeof(double) * nSize);

	return true;
}

// 将方阵初始化为单位矩阵
// 参数：
// 1. int nSize - 方阵行列数
// 返回值：bool 型，初始化是否成功
bool BDMatrix::MakeUnitMatrix(int nSize)
{
	if (!Init(nSize, nSize))
		return false;

	for (int i = 0; i < nSize; ++i)
		for (int j = 0; j < nSize; ++j)
			if (i == j)
				SetElement(i, j, 1);

	return true;
}

bool BDMatrix::clear()
{
	if (m_pData)
	{
		delete[] m_pData;
		m_pData = nullptr;
		m_nNumRows = 0;
		m_nNumColumns = 0;
	}
	return true;
}

// 设置矩阵各元素的值
// 参数：
// 1. double value[] - 一维数组，长度为m_nNumColumns*m_nNumRows，存储
// 矩阵各元素的值
// 返回值：无
void BDMatrix::SetData(double value[])
{
	// empty the memory;
	memset(m_pData, 0, sizeof(double) * m_nNumColumns * m_nNumRows);
	// copy data;
	memcpy(m_pData, value, sizeof(double) * m_nNumColumns * m_nNumRows);
}

// 设置指定元素的值
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
// 3. double value - 指定元素的值
// 返回值：bool 型，说明设置是否成功
bool BDMatrix::SetElement(int nRow, int nCol, double value)
{
	if (nCol < 0 || nCol >= m_nNumColumns || nRow < 0 || nRow >= m_nNumRows)
		return false;						// array bounds error;
	if (m_pData == nullptr)
		return false;							// bad pointer error;

	m_pData[nCol + nRow * m_nNumColumns] = value;

	return true;
}

bool BDMatrix::SetRow(int rowId, const std::vector<double>& row)
{
	for (size_t i = 0; i < row.size(); i++)
	{
		if (!SetElement(rowId, i, row[i]))
			return false;
	}
	return true;
}

// 设置指定元素的值
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
// 返回值：double 型，指定元素的值
double BDMatrix::GetElement(int nRow, int nCol) const
{
	assert(nCol >= 0 && nCol < m_nNumColumns && nRow >= 0 && nRow < m_nNumRows); // array bounds error;
	assert(m_pData);							// bad pointer error;
	return m_pData[nCol + nRow * m_nNumColumns];
}

// 获取矩阵的列数
// 参数：无
// 返回值：int 型，矩阵的列数
int	BDMatrix::GetNumColumns() const
{
	return m_nNumColumns;
}

// 获取矩阵的行数
// 参数：无
// 返回值：int 型，矩阵的行数
int	BDMatrix::GetNumRows() const
{
	return m_nNumRows;
}

// 获取矩阵的数据
// 参数：无
// 返回值：double型指针，指向矩阵各元素的数据缓冲区
double* BDMatrix::GetData() const
{
	return m_pData;
}

// 获取指定行的向量
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2.  double* pVector - 指向向量中各元素的缓冲区
// 返回值：int 型，向量中元素的个数，即矩阵的列数
int BDMatrix::GetRowVector(int nRow, double* pVector) const
{
	assert(pVector != nullptr);
	for (int j = 0; j < m_nNumColumns; ++j)
		pVector[j] = GetElement(nRow, j);

	return m_nNumColumns;
}

// 获取指定列的向量
// 参数：
// 1. int nCols - 指定的矩阵列数
// 2.  double* pVector - 指向向量中各元素的缓冲区
// 返回值：int 型，向量中元素的个数，即矩阵的行数
int BDMatrix::GetColVector(int nCol, double* pVector) const
{
	assert(pVector != nullptr);
	for (int i = 0; i < m_nNumRows; ++i)
		pVector[i] = GetElement(i, nCol);

	return m_nNumRows;
}

// 重载运算符=，给矩阵赋值
// 参数：
// 1. const BDMatrix& other - 用于给矩阵赋值的源矩阵
// 返回值：BDMatrix型的引用，所引用的矩阵与other相等
BDMatrix& BDMatrix::operator=(const BDMatrix& other)
{
	if (&other != this)
	{
		bool bSuccess = Init(other.GetNumRows(), other.GetNumColumns());
		assert(bSuccess);

		memcpy(m_pData, other.m_pData, sizeof(double) * m_nNumColumns * m_nNumRows);
	}
	return *this;
}

// 重载运算符==，判断矩阵是否相等
// 参数：
// 1. const BDMatrix& other - 用于比较的矩阵
// 返回值：bool 型，两个矩阵相等则为true，否则为false
bool BDMatrix::operator==(const BDMatrix& other) const
{
	// 首先检查行列数是否相等;
	if (m_nNumColumns != other.GetNumColumns() || m_nNumRows != other.GetNumRows())
		return false;

	for (int i = 0; i < m_nNumRows; ++i)
	{
		for (int j = 0; j < m_nNumColumns; ++j)
		{
			if (GetElement(i, j) != other.GetElement(i, j))
				return false;
		}
	}

	return true;
}

// 重载运算符!=，判断矩阵是否不相等
// 参数：
// 1. const BDMatrix& other - 用于比较的矩阵
// 返回值：bool 型，两个不矩阵相等则为true，否则为false
bool BDMatrix::operator!=(const BDMatrix& other) const
{
	return !(*this == other);
}

// 重载运算符+，实现矩阵的加法
// 参数：
// 1. const BDMatrix& other - 与指定矩阵相加的矩阵
// 返回值：BDMatrix型，指定矩阵与other相加之和
BDMatrix	BDMatrix::operator+(const BDMatrix& other) const
{
	// 首先检查行列数是否相等
	assert(m_nNumColumns == other.GetNumColumns() && m_nNumRows == other.GetNumRows());

	// 构造结果矩阵
	BDMatrix	result(*this);		// 拷贝构造
	// 矩阵加法
	for (int i = 0; i < m_nNumRows; ++i)
	{
		for (int j = 0; j < m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) + other.GetElement(i, j));
	}
	return result;
}

// 重载运算符-，实现矩阵的减法
// 参数：
// 1. const BDMatrix& other - 与指定矩阵相减的矩阵
// 返回值：BDMatrix型，指定矩阵与other相减之差
BDMatrix	BDMatrix::operator-(const BDMatrix& other) const
{
	// 首先检查行列数是否相等
	assert(m_nNumColumns == other.GetNumColumns() && m_nNumRows == other.GetNumRows());

	// 构造目标矩阵
	BDMatrix	result(*this);		// copy ourselves;
	// 进行减法操作
	for (int i = 0; i < m_nNumRows; ++i)
	{
		for (int j = 0; j < m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) - other.GetElement(i, j));
	}

	return result;
}

// 重载运算符*，实现矩阵的数乘
// 参数：
// 1. double value - 与指定矩阵相乘的实数
// 返回值：BDMatrix型，指定矩阵与value相乘之积
BDMatrix	BDMatrix::operator*(double value) const
{
	// 构造目标矩阵;
	BDMatrix	result(*this);		// copy ourselves
	// 进行数乘
	for (int i = 0; i < m_nNumRows; ++i)
	{
		for (int j = 0; j < m_nNumColumns; ++j)
			result.SetElement(i, j, result.GetElement(i, j) * value);
	}

	return result;
}

// 重载运算符*，实现矩阵的乘法
// 参数：
// 1. const BDMatrix& other - 与指定矩阵相乘的矩阵
// 返回值：BDMatrix型，指定矩阵与other相乘之积
BDMatrix	BDMatrix::operator*(const BDMatrix& other) const
{
	// 首先检查行列数是否符合要求
	assert(m_nNumColumns == other.GetNumRows());

	// construct the object we are going to return
	BDMatrix	result(m_nNumRows, other.GetNumColumns());

	// 矩阵乘法，即;
	//
	// [A][B][C]   [G][H]     [A*G + B*I + C*K][A*H + B*J + C*L];
	// [D][E][F] * [I][J] =   [D*G + E*I + F*K][D*H + E*J + F*L];
	//             [K][L];
	//
	double	value;
	for (int i = 0; i < result.GetNumRows(); ++i)
	{
		for (int j = 0; j < other.GetNumColumns(); ++j)
		{
			value = 0.0;
			for (int k = 0; k < m_nNumColumns; ++k)
			{
				value += GetElement(i, k) * other.GetElement(k, j);
			}

			result.SetElement(i, j, value);
		}
	}

	return result;
}

#ifdef inQT
void BDMatrix::load(const QString& path, bool existTitle)
{
	BDFile::importData(path, *this, existTitle);
}

void BDMatrix::save(QString path, std::vector<std::string> title)
{
	BDFile::exportData(path, *this, title);
}

void BDMatrix::print()
{
	forEachRow2([](int rowID, double* data, int ncol)
		{
			QString str;
			for (int i = 0; i < ncol - 1; i++)
			{
				str += QString::number(data[i]) + QString("  ");
			}
			str += QString::number(data[ncol - 1]);
			BDLOG str;
			//            BDLOGS rowID << str;
		});
}
#endif

double* BDMatrix::operator*(const double* other) const
{
	assert(m_nNumColumns == m_nNumRows);
	double* result = new double[m_nNumRows];
	double value;
	for (int i = 0; i < m_nNumRows; i++)
	{
		value = 0.0;
		for (int j = 0; j < m_nNumRows; j++)
			value += GetElement(i, j) * other[j];
		result[i] = value;
	}
	return result;
}

BDPoint3 BDMatrix::operator*(const BDPoint3 other) const
{
	assert(m_nNumColumns == m_nNumRows);
	assert(m_nNumColumns == 3);

	BDPoint3 result;
	double value;
	for (int i = 0; i < m_nNumColumns; i++)
	{
		value = 0.0;
		for (int j = 0; j < m_nNumColumns; j++)
			value += GetElement(i, j) * other[j];
		result[i] = value;
	}
	return result;
}

double* BDMatrix::operator*(const std::vector<double>& other) const
{
	assert(m_nNumColumns == m_nNumRows);
	double* result = new double[m_nNumRows];
	double value;
	for (int i = 0; i < m_nNumRows; i++)
	{
		value = 0.0;
		for (int j = 0; j < m_nNumRows; j++)
			value += GetElement(i, j) * other[j];
		result[i] = value;

	}
	return result;
}

// 重载运算符[],用于实现用[][]操作矩阵元素
double* BDMatrix::operator[](int nRow)
{
	assert(nRow >= 0 && nRow < m_nNumRows);
	//	if(nRow<0||nRow>=m_nNumRows);
	//		::AfxMessageBox("error");
	return m_pData + m_nNumColumns * nRow;
}

// 矩阵的转置
// 参数：无
// 返回值：BDMatrix型，指定矩阵转置矩阵
BDMatrix BDMatrix::Transpose() const
{
	BDMatrix m(m_nNumColumns, m_nNumRows);

	for (int i = 0; i < m_nNumRows; ++i)
		for (int j = 0; j < m_nNumColumns; ++j)
			m.SetElement(j, i, GetElement(i, j));

	return m;
}

// 实矩阵求逆的全选主元高斯－约当法
// 参数：无
// 返回值：bool型，求逆是否成功
bool BDMatrix::InvertGaussJordan()
{
	int* pnRow, * pnCol, i, j, k, l, u, v;
	double d = 0, p = 0;

	// 分配内存
	pnRow = new int[m_nNumColumns];
	pnCol = new int[m_nNumColumns];
	if (pnRow == nullptr || pnCol == nullptr)
		return false;

	// 消元
	for (k = 0; k <= m_nNumColumns - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= m_nNumColumns - 1; i++)
		{
			for (j = k; j <= m_nNumColumns - 1; j++)
			{
				l = i * m_nNumColumns + j; p = fabs(m_pData[l]);
				if (p > d)
				{
					d = p;
					pnRow[k] = i;
					pnCol[k] = j;
				}
			}
		}

		// 失败
		if (d == 0.0)
		{
			delete[] pnRow;
			delete[] pnCol;
			return false;
		}

		if (pnRow[k] != k)
		{
			for (j = 0; j <= m_nNumColumns - 1; j++)
			{
				u = k * m_nNumColumns + j;
				v = pnRow[k] * m_nNumColumns + j;
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		if (pnCol[k] != k)
		{
			for (i = 0; i <= m_nNumColumns - 1; i++)
			{
				u = i * m_nNumColumns + k;
				v = i * m_nNumColumns + pnCol[k];
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		l = k * m_nNumColumns + k;
		m_pData[l] = 1.0 / m_pData[l];
		for (j = 0; j <= m_nNumColumns - 1; j++)
		{
			if (j != k)
			{
				u = k * m_nNumColumns + j;
				m_pData[u] = m_pData[u] * m_pData[l];
			}
		}

		for (i = 0; i <= m_nNumColumns - 1; i++)
		{
			if (i != k)
			{
				for (j = 0; j <= m_nNumColumns - 1; j++)
				{
					if (j != k)
					{
						u = i * m_nNumColumns + j;
						m_pData[u] = m_pData[u] - m_pData[i * m_nNumColumns + k] * m_pData[k * m_nNumColumns + j];
					}
				}
			}
		}

		for (i = 0; i <= m_nNumColumns - 1; i++)
		{
			if (i != k)
			{
				u = i * m_nNumColumns + k;
				m_pData[u] = -m_pData[u] * m_pData[l];
			}
		}
	}

	// 调整恢复行列次序
	for (k = m_nNumColumns - 1; k >= 0; k--)
	{
		if (pnCol[k] != k)
		{
			for (j = 0; j <= m_nNumColumns - 1; j++)
			{
				u = k * m_nNumColumns + j;
				v = pnCol[k] * m_nNumColumns + j;
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}

		if (pnRow[k] != k)
		{
			for (i = 0; i <= m_nNumColumns - 1; i++)
			{
				u = i * m_nNumColumns + k;
				v = i * m_nNumColumns + pnRow[k];
				p = m_pData[u];
				m_pData[u] = m_pData[v];
				m_pData[v] = p;
			}
		}
	}

	// 清理内存
	delete[] pnRow;
	delete[] pnCol;

	// 成功返回
	return true;
}

BDMatrix& BDMatrix::combineH(BDMatrix& other)
{
	assert(GetNumRows() == other.GetNumRows());

	int nrow = this->GetNumRows();
	int ncol = this->GetNumColumns() + other.GetNumColumns();

	BDMatrix m(nrow, ncol);
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < ncol; j++)
		{
			double v;
			if (j < GetNumColumns())
				v = GetElement(i, j);
			else
				v = other.GetElement(i, j - GetNumColumns());
			m.SetElement(i, j, v);
		}
	}
	*this = m;
	return *this;
}

BDMatrix& BDMatrix::combineV(BDMatrix& other)
{
	assert(GetNumColumns() == other.GetNumColumns());

	int nrow = this->GetNumRows() + other.GetNumRows();
	int ncol = this->GetNumColumns();

	BDMatrix m(nrow, ncol);
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < ncol; j++)
		{
			double v;
			if (i < GetNumRows())
				v = GetElement(i, j);
			else
				v = other.GetElement(i - GetNumRows(), j);
			m.SetElement(i, j, v);
		}
	}
	*this = m;
	return *this;
}

void BDMatrix::trans(double* source, double* desm) //矩阵转换 desm = M*source
{
	assert(m_nNumColumns == m_nNumRows);

	double value;
	for (int i = 0; i < m_nNumRows; i++)
	{
		value = 0.0;
		for (int j = 0; j < m_nNumRows; j++)
			value += GetElement(i, j) * source[j];
		desm[i] = value;
	}
}

void BDMatrix::toVector(std::vector<std::vector<double> >& value)
{
	value.clear();
	std::vector<double> row;
	row.resize(GetNumColumns());
	value.resize(GetNumRows(), row);
	forEachRow2([&value](int rowId, double* row, int len)
		{
			for (int c = 0; c < len; c++)
			{
				value[rowId][c] = row[c];
			}
		});
}

void BDMatrix::forEachElement(std::function<void(int, int, double)> op)
{
	for (int i = 0; i < GetNumRows(); i++)
	{
		for (int j = 0; j < GetNumColumns(); j++)
		{
			op(i, j, (*this)[i][j]);
		}
	}
}

void BDMatrix::forEachRow(std::function<void(double*, int)> op)
{
	for (int i = 0; i < GetNumRows(); i++)
	{
		op((*this)[i], GetNumColumns());
	}
}

void BDMatrix::forEachRow2(std::function<void(int, double*, int)> op)
{
	for (int i = 0; i < GetNumRows(); i++)
	{
		op(i, (*this)[i], GetNumColumns());
	}
}

BDMatrix BDMatrix::colume(const BDVector& cols)
{
	BDMatrix p(GetNumRows(), cols.size());
	forEachRow2([&](int rowId, double* data, int ncol)
		{
			int id = 0;
			for (int colId : cols) {
				p[rowId][id] = data[colId];
				id++;
			}
		});
	return p;
}

BDVector BDMatrix::colume(unsigned col)
{
	assert(col < GetNumColumns());
	BDVector v;
	v.resize(GetNumRows());
	for (int i = 0; i < GetNumRows(); i++) {
		v[i] = (*this)[i][col];
	}
	return v;
}

double* BDMatrix::firstRow()
{
	return (*this)[0];
}

double* BDMatrix::lastRow()
{
	return (*this)[GetNumRows() - 1];
}

/****************************************************
* BDString
* *************************************************/

BDString::BDString()
{}

void BDString::replaceAll(const std::string& o, const std::string& n)
{
	for (string::size_type pos(0); pos != string::npos; pos += n.length())
	{
		pos = this->find(o, pos);
		if (pos == string::npos)
			break;
		this->replace(pos, o.length(), n);
	}
}

/***********************************************************************
BDCoordTransTool
***********************************************************************/

BDCoordTransTool::BDCoordTransTool()
{}

BDCoordTransTool::~BDCoordTransTool()
{
}

BDPoint3 BDCoordTransTool::LBH2Earth(BDPoint3 lla)
{
	// 《总体设计(上)》Eq(4.10-7)，p411
	double lng = lla[0];
	double lat = lla[1];
	double h = lla[2];

	if (lng > 180 / a57) lng = 180 / a57;
	if (lng < -180 / a57) lng = -180 / a57;

	if (lat > 90 / a57) lat = 90 / a57;
	if (lat < -90 / a57) lat = -90 / a57;

	double a = BDEarth::a, b = BDEarth::b;
	if (!IS_ELLIPSOID) {
		a = BDEarth::ro;
		b = BDEarth::ro;
	}

	double e2 = (pow(a, 2.0) - pow(b, 2.0)) / pow(a, 2.0);
	double w = pow(1 - e2 * pow(sin(lat), 2.0), 0.5);
	double N = a / w;

	BDPoint3 r;
	r.x = (N + h) * cos(lat) * cos(lng);
	r.y = (N + h) * cos(lat) * sin(lng);
	r.z = (N * (1 - e2) + h) * sin(lat);
	return r;
}

BDPoint3 BDCoordTransTool::LBH2Earthd(BDPoint3 lla)
{
	lla.deg2rad();
	return BDCoordTransTool::LBH2Earth(lla);
}

void BDCoordTransTool::LBH2Earth(BDMatrix& lla)
{
	lla.forEachRow([](double* row, int ncol)
		{
			BDPoint3 p(row[0], row[1], row[2]);
			p = BDCoordTransTool::LBH2Earth(p);
			row[0] = p.x;
			row[1] = p.y;
			row[2] = p.z;
		});
}

void BDCoordTransTool::LBH2Earthd(BDMatrix& lla)
{
	lla.forEachRow([](double* row, int ncol)
		{
			BDPoint3 p(row[0], row[1], row[2]);
			p = BDCoordTransTool::LBH2Earthd(p);
			row[0] = p.x;
			row[1] = p.y;
			row[2] = p.z;
		});
}

double BDCoordTransTool::B2Phi(double B)
{
	return atan((1 - BDEarth::e2) * tan(B));
}

double BDCoordTransTool::Phi2B(double phi)
{
	return atan(tan(phi) / (1 - BDEarth::e2));
}

/**
* @brief 参考飞行动力学p119 弹道导弹弹道学p17
* @param x0
* @param ellipsoid
* @return
*/
BDPoint3 BDCoordTransTool::Earth2LPH(BDPoint3 x0)
{
	// 地心地固系坐标
	double x = x0.x;
	double y = x0.y;
	double z = x0.z;

	// 地心距
	double r = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));

	// 经度
	double lambda = iszero(y) ? 0 : y / fabs(y) * (PI / 2 - asin(x / sqrt(pow(x, 2.0) + pow(y, 2.0))));

	// 地心纬度
	double phi = asin(z / r);

	// 高度
	double h = 0.0;
	if (IS_ELLIPSOID) {
		h = r - BDEarth::calR0(phi);
	}
	else
		h = r - BDEarth::ro;
	return BDPoint3(lambda, phi, h);
}

BDPoint3 BDCoordTransTool::Earth2LPHd(BDPoint3 x0)
{
	auto x = BDCoordTransTool::Earth2LPH(x0);
	x.rad2deg();
	return x;
}

void BDCoordTransTool::Earth2LPH(BDMatrix& x0)
{
	x0.forEachRow([](double* row, int ncol)
		{
			BDPoint3 p(row[0], row[1], row[2]);
			p = BDCoordTransTool::Earth2LPH(p);
			row[0] = p.x;
			row[1] = p.y;
			row[2] = p.z;
		});
}

void BDCoordTransTool::Earth2LPHd(BDMatrix& x0)
{
	x0.forEachRow([](double* row, int ncol)
		{
			BDPoint3 p(row[0], row[1], row[2]);
			p = BDCoordTransTool::Earth2LPHd(p);
			row[0] = p.x;
			row[1] = p.y;
			row[2] = p.z;
		});
}

double BDCoordTransTool::calDistance(BDPoint3 LLA1, BDPoint3 LLA2)
{
	LLA1 = LBH2Earth(LLA1);
	LLA2 = LBH2Earth(LLA2);
	return LLA1.Distance(LLA2);
}

BDVector BDCoordTransTool::calTwoPointRelation(double lambda1, double B1, double lambda2, double B2)
{
	// 计算地心纬度
	double phi1, phi2;
	if (IS_ELLIPSOID)
	{
		phi1 = atan((1 - BDEarth::e2) * tan(B1));
		phi2 = atan((1 - BDEarth::e2) * tan(B2));
	}
	else
	{
		phi1 = B1;
		phi2 = B2;
	}

	//由东至西跨越正负180度经度线的判据（跨越线的经度差小于180°）
	bool EastToWestCross180Lon = (fabs(lambda2 - lambda1 + 2 * PI) < PI);

	//由西至东跨越正负180度经度线的判据（跨越线的经度差小于180°）
	bool WestToEastCross180Lon = (fabs(lambda2 - lambda1 - 2 * PI) < PI);

	//计算经度差
	double dlon;
	if (EastToWestCross180Lon)
	{
		dlon = lambda2 - lambda1 + 2 * PI;
	}
	else if (WestToEastCross180Lon)
	{
		dlon = lambda2 - lambda1 - 2 * PI;
	}
	else
	{
		dlon = lambda2 - lambda1;
	}

	//计算(lon1,lat1)到(lon2,lat2)的射程角
	double PHIs = acos(cos(phi1) * cos(phi2) * cos(dlon) + sin(phi1) * sin(phi2));

	double angleRange; // 角射程
	double lineRange;  // 线射程
	double azimuth;    // 球面方位角
	double PRECISION = 1e-10;
	if (fabs(PHIs) >= PRECISION * 1e5)
	{
		//定义发射点目标点球面方位角余弦值
		double COSPSIs = (sin(phi2) - sin(phi1) * cos(PHIs)) / (cos(phi1) * sin(PHIs));

		//定义发射点目标点球面方位角正弦值
		double SINPSIs = cos(phi2) * sin(dlon) / sin(PHIs);

		//定义发射点目标点球面方位角
		double PSIs;
		if (fabs(COSPSIs) > (1 - PRECISION * 10e2))
		{
			COSPSIs = BDMath::sign(COSPSIs);
		}

		//
		if (SINPSIs >= 0)
			PSIs = acos(COSPSIs);
		else
			PSIs = -acos(COSPSIs);

		if (PSIs < 0)
			PSIs = PSIs + 2 * PI;

		//计算两点角射程
		angleRange = PHIs;

		//计算两点线射程
		lineRange = BDEarth::ro * PHIs;

		//计算(lon1,lat1)到(lon2,lat2)的球面方位角
		azimuth = PSIs;
	}
	else
	{
		//计算两点角射程
		angleRange = PHIs;

		//计算两点线射程
		lineRange = BDEarth::ro * PHIs;

		//计算(lon1,lat1)到(lon2,lat2)的球面方位角
		azimuth = 0;
	}

	//返回两点球面关系结构
	return BDVector({ angleRange, lineRange, azimuth });
}

double BDCoordTransTool::calS(double lambda1, double B1, double lambda2, double B2)
{
	return calTwoPointRelation(lambda1, B1, lambda2, B2)[1];
}

BDPoint3 BDCoordTransTool::Earth2Launch(const BDVector& launchInfo, BDPoint3 p)
{
	double lambda0 = launchInfo[0];
	double phi0 = launchInfo[1];
	double h0 = launchInfo[2];
	double alpha0 = launchInfo[3];

	BDPoint3 lla(lambda0, phi0, h0);
	BDPoint3 x0 = LBH2Earth(lla);

	BDMatrix m(3, 3);
	m.SetElement(0, 0, -sin(alpha0) * sin(lambda0) - cos(alpha0) * sin(phi0) * cos(lambda0));
	m.SetElement(1, 0, cos(phi0) * cos(lambda0));
	m.SetElement(2, 0, -cos(alpha0) * sin(lambda0) + sin(alpha0) * sin(phi0) * cos(lambda0));
	m.SetElement(0, 1, sin(alpha0) * cos(lambda0) - cos(alpha0) * sin(phi0) * sin(lambda0));
	m.SetElement(1, 1, cos(phi0) * sin(lambda0));
	m.SetElement(2, 1, cos(alpha0) * cos(lambda0) + sin(alpha0) * sin(phi0) * sin(lambda0));
	m.SetElement(0, 2, cos(alpha0) * cos(phi0));
	m.SetElement(1, 2, sin(phi0));
	m.SetElement(2, 2, -sin(alpha0) * cos(phi0));

	return m * (p - x0);
}

BDPoint3 BDCoordTransTool::Earth2Launchv(const BDVector& launchInfo, BDPoint3 p)
{
	double lambda0 = launchInfo[0];
	double phi0 = launchInfo[1];
	double h0 = launchInfo[2];
	double alpha0 = launchInfo[3];

	BDMatrix m(3, 3);
	m.SetElement(0, 0, -sin(alpha0) * sin(lambda0) - cos(alpha0) * sin(phi0) * cos(lambda0));
	m.SetElement(1, 0, cos(phi0) * cos(lambda0));
	m.SetElement(2, 0, -cos(alpha0) * sin(lambda0) + sin(alpha0) * sin(phi0) * cos(lambda0));
	m.SetElement(0, 1, sin(alpha0) * cos(lambda0) - cos(alpha0) * sin(phi0) * sin(lambda0));
	m.SetElement(1, 1, cos(phi0) * sin(lambda0));
	m.SetElement(2, 1, cos(alpha0) * cos(lambda0) + sin(alpha0) * sin(phi0) * sin(lambda0));
	m.SetElement(0, 2, cos(alpha0) * cos(phi0));
	m.SetElement(1, 2, sin(phi0));
	m.SetElement(2, 2, -sin(alpha0) * cos(phi0));

	return m * p;
}

BDPoint3 BDCoordTransTool::Launch2Earth(const BDVector& launchInfo, BDPoint3 p)
{
	double lambda0 = launchInfo[0];
	double B0 = launchInfo[1];
	double h0 = launchInfo[2];
	double A0 = launchInfo[3];

	BDPoint3 LBH(lambda0, B0, h0);
	BDPoint3 x0 = LBH2Earth(LBH);

	BDMatrix m(3, 3);
	m.SetElement(0, 0, -sin(A0) * sin(lambda0) - cos(A0) * sin(B0) * cos(lambda0));
	m.SetElement(1, 0, cos(B0) * cos(lambda0));
	m.SetElement(2, 0, -cos(A0) * sin(lambda0) + sin(A0) * sin(B0) * cos(lambda0));
	m.SetElement(0, 1, sin(A0) * cos(lambda0) - cos(A0) * sin(B0) * sin(lambda0));
	m.SetElement(1, 1, cos(B0) * sin(lambda0));
	m.SetElement(2, 1, cos(A0) * cos(lambda0) + sin(A0) * sin(B0) * sin(lambda0));
	m.SetElement(0, 2, cos(A0) * cos(B0));
	m.SetElement(1, 2, sin(B0));
	m.SetElement(2, 2, -sin(A0) * cos(B0));
	m = m.Transpose();

	return m * p + x0;
}

void BDCoordTransTool::Launch2ECI(const BDVector& launchInfo, double* in)
{
	auto ecef = BDCoordTransTool::Launch2Earth(launchInfo, BDPoint3(in[1], in[2], in[3]));
	auto eci = BDCoordTransTool::Earth2EarthI(in[0], ecef);
	auto ecefv = BDCoordTransTool::Launch2Earthv(launchInfo, BDPoint3(in[4], in[5], in[6]));
	auto eciv = BDCoordTransTool::Earth2EarthIv(in[0], ecef, ecefv);
	in[1] = eci[0];
	in[2] = eci[1];
	in[3] = eci[2];
	in[4] = eciv[0];
	in[5] = eciv[1];
	in[6] = eciv[2];
}

void BDCoordTransTool::Launch2ECI(const BDVector& launchInfo, BDVector& in)
{
	auto ecef = BDCoordTransTool::Launch2Earth(launchInfo, BDPoint3(in[1], in[2], in[3]));
	auto eci = BDCoordTransTool::Earth2EarthI(in[0], ecef);
	auto ecefv = BDCoordTransTool::Launch2Earthv(launchInfo, BDPoint3(in[4], in[5], in[6]));
	auto eciv = BDCoordTransTool::Earth2EarthIv(in[0], ecef, ecefv);
	in[1] = eci[0];
	in[2] = eci[1];
	in[3] = eci[2];
	in[4] = eciv[0];
	in[5] = eciv[1];
	in[6] = eciv[2];
}

//    BDPoint3 BDCoordTransTool::Launch2ECI(const BDVector &launchInfo, double t, const BDPoint3 &p, const BDPoint3 &v, BDVector &result)
//    {
//        auto ecef = BDCoordTransTool::Launch2Earth(launchInfo, p);
//        auto eci = BDCoordTransTool::Earth2EarthI(t, ecef);

//        auto ecefv = BDCoordTransTool::Launch2Earthv(launchInfo, v);
//        auto eciv = BDCoordTransTool::Earth2EarthIv(t, ecef, ecefv);

//        result.resize(7);
//        result[0] = t;

//        result[1] = eci[0];
//        result[2] = eci[1];
//        result[3] = eci[2];

//        result[4] = eciv[0];
//        result[5] = eciv[1];
//        result[6] = eciv[2];
//    }


BDPoint3 BDCoordTransTool::Launch2LPH(const BDVector& launchInfo, const BDPoint3& p)
{
	auto ecef = BDCoordTransTool::Launch2Earth(launchInfo, p);
	return BDCoordTransTool::Earth2LPH(ecef);
}

BDPoint3 BDCoordTransTool::Launch2Earthv(const BDVector& launchInfo, BDPoint3 v)
{
	double lambda0 = launchInfo[0];
	double B0 = launchInfo[1];
	double A0 = launchInfo[3];

	BDMatrix m(3, 3);
	m.SetElement(0, 0, -sin(A0) * sin(lambda0) - cos(A0) * sin(B0) * cos(lambda0));
	m.SetElement(1, 0, cos(B0) * cos(lambda0));
	m.SetElement(2, 0, -cos(A0) * sin(lambda0) + sin(A0) * sin(B0) * cos(lambda0));
	m.SetElement(0, 1, sin(A0) * cos(lambda0) - cos(A0) * sin(B0) * sin(lambda0));
	m.SetElement(1, 1, cos(B0) * sin(lambda0));
	m.SetElement(2, 1, cos(A0) * cos(lambda0) + sin(A0) * sin(B0) * sin(lambda0));
	m.SetElement(0, 2, cos(A0) * cos(B0));
	m.SetElement(1, 2, sin(B0));
	m.SetElement(2, 2, -sin(A0) * cos(B0));

	m = m.Transpose();
	return m * v;
}

BDPoint3 BDCoordTransTool::Body2Launchv(const BDPoint3& phi_psi_gamma, const BDPoint3& p)
{
	double phi = phi_psi_gamma[0];
	double psi = phi_psi_gamma[1];
	double gamma = phi_psi_gamma[2];

	BDMatrix m(3, 3);
	m.SetElement(0, 0, cos(phi) * cos(psi));
	m.SetElement(0, 1, cos(phi) * sin(psi) * sin(gamma) - sin(phi) * cos(gamma));
	m.SetElement(0, 2, cos(phi) * sin(psi) * cos(gamma) + sin(phi) * sin(gamma));

	m.SetElement(1, 0, sin(phi) * cos(psi));
	m.SetElement(1, 1, sin(phi) * sin(psi) * sin(gamma) + cos(phi) * cos(gamma));
	m.SetElement(1, 2, sin(phi) * sin(psi) * cos(gamma) - cos(phi) * sin(gamma));

	m.SetElement(2, 0, -sin(psi));
	m.SetElement(2, 1, cos(psi) * sin(gamma));
	m.SetElement(2, 2, cos(psi) * cos(gamma));

	return m * p;
}

BDPoint3 BDCoordTransTool::Body2Launchp(const BDPoint3& body, const BDPoint3& phi_psi_gamma, const BDPoint3& p)
{
	return Body2Launchv(phi_psi_gamma, p) + body;
}

BDPoint3 BDCoordTransTool::Launch2Body(BDPoint3 phi_psi_gamma, BDPoint3 p, BDPoint3 body)
{
	double phi = phi_psi_gamma[0];
	double psi = phi_psi_gamma[1];
	double gamma = phi_psi_gamma[2];

	BDMatrix m(3, 3);
	m.SetElement(0, 0, cos(phi) * cos(psi));
	m.SetElement(0, 1, cos(phi) * sin(psi) * sin(gamma) - sin(phi) * cos(gamma));
	m.SetElement(0, 2, cos(phi) * sin(psi) * cos(gamma) + sin(phi) * sin(gamma));
	m.SetElement(1, 0, sin(phi) * cos(psi));
	m.SetElement(1, 1, sin(phi) * sin(psi) * sin(gamma) + cos(phi) * cos(gamma));
	m.SetElement(1, 2, sin(phi) * sin(psi) * cos(gamma) - cos(phi) * sin(gamma));
	m.SetElement(2, 0, -sin(psi));
	m.SetElement(2, 1, cos(psi) * sin(gamma));
	m.SetElement(2, 2, cos(psi) * cos(gamma));
	m = m.Transpose();

	return m * (p - body);
}

BDPoint3 BDCoordTransTool::Earth2ENU(BDPoint3 ecef, BDPoint3 radarLBH)
{
	auto radarECEF = BDCoordTransTool::LBH2Earth(radarLBH);
	double L = radarLBH[0]; // 经度
	double B = radarLBH[1]; // 纬度

	BDMatrix m(3, 3);
	m.SetRow(0, { -sin(L), cos(L), 0 });
	m.SetRow(1, { -cos(L) * sin(B), -sin(L) * sin(B), cos(B) });
	m.SetRow(2, { cos(L) * cos(B), sin(L) * cos(B), sin(B) });

	return m * (ecef - radarECEF);
}

BDPoint3 BDCoordTransTool::Earth2EarthI(double t, const BDPoint3& x)
{
	double a = BDEarth::OMEAGA * t;
	BDMatrix m(3, 3);
	m[0][0] = cos(a);     m[0][1] = -sin(a);      m[0][2] = 0.0;
	m[1][0] = sin(a);     m[1][1] = cos(a);       m[1][2] = 0.0;
	m[2][0] = 0.0;        m[2][1] = 0.0;          m[2][2] = 1.0;
	return m * x;
}

BDPoint3 BDCoordTransTool::Earth2EarthIv(double t, const BDPoint3& x, const BDPoint3& v)
{
	double a = BDEarth::OMEAGA * t;
	BDMatrix A(3, 3);
	A[0][0] = cos(a);     A[0][1] = -sin(a);      A[0][2] = 0.0;
	A[1][0] = sin(a);     A[1][1] = cos(a);       A[1][2] = 0.0;
	A[2][0] = 0.0;        A[2][1] = 0.0;          A[2][2] = 1.0;

	BDMatrix B(3, 3);
	B[0][0] = -sin(a);    B[0][1] = -cos(a);     B[0][2] = 0.0;
	B[1][0] = cos(a);     B[1][1] = -sin(a);     B[1][2] = 0.0;
	B[2][0] = 0.0;        B[2][1] = 0.0;         B[2][2] = 0.0;
	B = B * BDEarth::OMEAGA;

	return A * v + B * x;
}

BDPoint3 BDCoordTransTool::EarthI2Earth(double t, const BDPoint3& x)
{
	double a = BDEarth::OMEAGA * t;
	BDMatrix m(3, 3);
	m[0][0] = cos(a);
	m[0][1] = sin(a);
	m[0][2] = 0.0;
	m[1][0] = -sin(a);
	m[1][1] = cos(a);
	m[1][2] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;
	return m * x;
}


BDPoint3 BDCoordTransTool::Velocity2Body(const BDPoint3& x, const BDPoint3& alpha_beta)
{
	double alpha = alpha_beta[0];
	double beta = alpha_beta[1];
	BDMatrix m(3, 3);
	m.SetRow(0, { cos(alpha) * cos(beta),    sin(alpha),    -cos(alpha) * sin(beta) });
	m.SetRow(1, { -sin(alpha) * cos(beta),   cos(alpha),    sin(alpha) * sin(beta) });
	m.SetRow(2, { sin(beta),                 0,             cos(beta) });
	return m * x;
}

BDPoint3 BDCoordTransTool::Velocity2Launch(const BDPoint3& x, const BDPoint3& theta_sigma_gammav)
{
	double theta = theta_sigma_gammav[0];
	double sigma = theta_sigma_gammav[1];
	double gammav = theta_sigma_gammav[2];

	BDMatrix m(3, 3);
	m[0][0] = cos(theta) * cos(sigma);
	m[1][0] = sin(theta) * cos(sigma);
	m[2][0] = -sin(sigma);

	m[0][1] = cos(theta) * sin(sigma) * sin(gammav) - sin(theta) * cos(gammav);
	m[1][1] = sin(theta) * sin(sigma) * sin(gammav) + cos(theta) * cos(gammav);
	m[2][1] = cos(sigma) * sin(gammav);

	m[0][2] = cos(theta) * sin(sigma) * cos(gammav) + sin(theta) * sin(gammav);
	m[1][2] = sin(theta) * sin(sigma) * cos(gammav) - cos(theta) * sin(gammav);
	m[2][2] = cos(sigma) * cos(gammav);

	return m * x;
}

BDPoint3 BDCoordTransTool::Launch2LaunchIv(double B0, double A0, double t, const BDPoint3& v)
{
	double X = cos(B0) * cos(A0);
	double Y = sin(B0);
	double Z = -cos(B0) * sin(A0);
	double S = sin(BDEarth::OMEAGA * t);
	double C = cos(BDEarth::OMEAGA * t);

	BDMatrix m1(3, 3);   // 发射到发惯系转换矩阵
	m1[0][0] = X * X * (1 - C) + C;
	m1[0][1] = X * Y * (1 - C) - Z * S;
	m1[0][2] = X * Z * (1 - C) + Y * S;

	m1[1][0] = X * Y * (1 - C) + Z * S;
	m1[1][1] = Y * Y * (1 - C) + C;
	m1[1][2] = Y * Z * (1 - C) - X * S;

	m1[2][0] = X * Z * (1 - C) - Y * S;
	m1[2][1] = Y * Z * (1 - C) + X * S;
	m1[2][2] = Z * Z * (1 - C) + C;

	return m1 * v;
}

BDPoint3 BDCoordTransTool::Launch2LaunchIp(double B0, double A0, double t,
	const BDPoint3& p, const BDPoint3& R0)
{
	return Launch2LaunchIv(B0, A0, t, p + R0) - R0;
}

BDVector BDCoordTransTool::Launch2LaunchIpv(double B0, double A0, double t,
	const BDPoint3& p, const BDPoint3& v,
	const BDPoint3& w, const BDPoint3& R0)
{
	auto ra = Launch2LaunchIv(B0, A0, t, p + R0);
	auto v1 = Launch2LaunchIv(B0, A0, t, v);
	auto va = v1 + w.cross(ra);
	auto pa = ra - R0;
	return BDVector({ pa.x, pa.y, pa.z, va.x, va.y, va.z });
}

BDPoint3 BDCoordTransTool::movePoint(const BDPoint3& centerLLA, double alpha, double xDis, double zDis)
{
	double x = -zDis * sin(alpha) + xDis * cos(alpha);
	double z = -zDis * cos(alpha) - xDis * sin(alpha);
	double y = 0;

	double a = BDEarth::a;
	double b = BDEarth::b;
	double e2 = BDEarth::e2;

	double Fio = atan((1.0 - e2) * tan(centerLLA.y));
	double Ro = centerLLA.z + a * b / sqrt(a * a * sin(Fio) * sin(Fio) + b * b * cos(Fio) * cos(Fio));
	double Uo = centerLLA.y - Fio;
	double Rox = -Ro * sin(Uo) * cos(0.0);
	double Roy = Ro * cos(Uo);
	double Roz = Ro * sin(Uo) * sin(0.0);

	double r = sqrt((x + Rox) * (x + Rox) + (y + Roy) * (y + Roy) + (z + Roz) * (z + Roz));
	double Fi = asin(((x + Rox) * cos(0.0) * cos(centerLLA.y) + (y + Roy) * sin(centerLLA.y) - (z + Roz) * sin(0.0) * cos(centerLLA.y)) / r);
	double Hdxd = r - b / sqrt(1.0 - e2 * cos(Fi) * cos(Fi));
	double Bdxd = atan(tan(Fi) / (1 - e2));
	double yd = -sin(centerLLA.y) * cos(0.0) * (x + Rox) + cos(centerLLA.y) * (y + Roy) + sin(centerLLA.y) * sin(0.0) * (z + Roz);
	double zd = sin(0.0) * (x + Rox) + cos(0.0) * (z + Roz);

	double Ldxd;
	if (yd > 0)
		Ldxd = centerLLA.x + atan(zd / yd);
	else
		Ldxd = centerLLA.x + PI + atan(zd / yd);

	BDPoint3 newPoint(Ldxd, Bdxd, centerLLA.z);
	return newPoint;
}


/***********************************************************************
BDFileReader
***********************************************************************/

BDFile::BDFile()
{
}

BDFile::~BDFile()
{
}

#ifdef inQT
bool BDFile::importData(const QString filePath, BDMatrix& data, bool existTitle)
{
	QFile resultFile(filePath);

	if (!resultFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
		BDLOG QStringLiteral("错误：文件打开失败。") << filePath;
		return false;
	}

	// 表头
	if (existTitle)
		resultFile.readLine();

	QChar splitChar(' ');
	int ncol = 0;
	int nrow = 0;
	while (!resultFile.atEnd())
	{
		QStringList items = QString(resultFile.readLine()).replace('\t', splitChar).replace(',', splitChar).trimmed().split(splitChar);
		items.removeAll("");

		if (ncol == 0)
			ncol = items.size();

		if (items.size() != ncol)
			continue;
		nrow++;
	}
	data.Init(nrow, ncol);
	resultFile.close();

	// 读取数据内容
	resultFile.open(QIODevice::ReadOnly | QIODevice::Text);
	if (existTitle)
		resultFile.readLine();

	int r = 0;
	while (!resultFile.atEnd())
	{
		QStringList items = QString(resultFile.readLine()).replace('\t', splitChar).replace(',', splitChar).trimmed().split(splitChar);
		items.removeAll("");
		if (items.size() != ncol)
			continue;

		for (int i = 0; i < ncol; i++)
		{
			data.SetElement(r, i, items[i].toDouble());
		}
		r++;
	}
	resultFile.close();
	return true;
}

bool BDFile::exportData(const QString filePath, BDMatrix& data, const std::vector<string>& title)
{
	// 创建新文件
	QFile file(filePath);
	file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate);
	QTextStream out(&file);
	out.setRealNumberPrecision(10);

	char split = ',';

	if (title.size() > 0)
	{
		for (int i = 0; i < title.size() - 1; i++)
			out << QString::fromStdString(title[i]) << split;
		out << QString::fromStdString(title[title.size() - 1]) << "\n";
	}

	for (int i = 0; i < data.GetNumRows(); i++)
	{
		for (int j = 0; j < data.GetNumColumns(); j++)
		{
			out << data[i][j];
			if (j < data.GetNumColumns() - 1)
				out << split;
		}
		if (i < data.GetNumRows() - 1)
			out << "\n";
	}
	file.close();
	return true;
}

bool BDFile::appendData(const QString filePath, BDMatrix& data, const QStringList& title)
{
	// 创建新文件
	QFile file(filePath);
	file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append);
	QTextStream out(&file);
	out.setRealNumberPrecision(10);
	out << "\n";

	if (title.size() > 0)
	{
		for (int i = 0; i < title.size() - 1; i++)
			out << title[i] << " ";
		out << title[title.size() - 1] << "\n";
	}

	for (int i = 0; i < data.GetNumRows(); i++)
	{
		for (int j = 0; j < data.GetNumColumns(); j++)
		{
			out << data[i][j];
			if (j < data.GetNumColumns() - 1)
				out << " ";
		}
		if (i < data.GetNumRows() - 1)
			out << "\n";
	}
	file.close();
	return true;
}

int BDFile::countLines(QString filePath)
{
	QFile resultFile(filePath);

	int c = 0;
	if (!resultFile.open(QIODevice::ReadOnly | QIODevice::Text))
		return c;

	while (!resultFile.atEnd())
	{
		resultFile.readLine();
		c++;
	}
	return c;
}

bool BDFile::toFile(QString filePath, QString content)
{
	ofstream outFile;
	outFile.open(filePath.toStdString(), ios::out);
	if (outFile) {
		outFile << content.toStdString();
		outFile.close();
		return true;
	}
	else
		return false;
	//        if(!BDFileReader::removeFile(filePath))
	//            return false;

	//        QFile file(filePath);
	//        if(!file.open(QIODevice::WriteOnly|QIODevice::Text|QIODevice::Truncate))
	//        {
	//            return false;
	//        }
	//        QTextStream out(&file);
	//        out<<content.toUtf8() <<endl;
	//        file.close();
	//        return true;
}

bool BDFile::readFile(QString filePath, QString& content)
{
	QFile file(filePath);

	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return false;

	content.clear();

	while (!file.atEnd()) {
		content = QString(file.readAll());
		QByteArray bytes = content.toUtf8();  //获取字节
		content = QString(bytes);
	}
	return true;
}

bool BDFile::removeDir(QString dirName)
{
	QDir dir(dirName);

	QString tmpdir = "";
	if (!dir.exists()) {
		return false;
	}

	QFileInfoList fileInfoList = dir.entryInfoList();
	foreach(QFileInfo fileInfo, fileInfoList) {
		if (fileInfo.fileName() == "." || fileInfo.fileName() == "..")
			continue;

		if (fileInfo.isDir()) {
			tmpdir = dirName + ("/") + fileInfo.fileName();
			removeDir(tmpdir);    //递归删除
			dir.rmdir(fileInfo.fileName()); //移除子目录
		}
		else if (fileInfo.isFile()) {
			QFile tmpFile(fileInfo.fileName());
			tmpFile.setPermissions(QFile::WriteOwner);
			dir.remove(tmpFile.fileName()); //删除文件
		}
	}

	dir.cdUp(); //返回上级目录，因为只有返回上级目录，才可以删除这个目录
	if (dir.exists(dirName)) {
		if (!dir.rmdir(dirName))
			return false;
	}
	return true;

}

bool BDFile::removeFile(QString filename)
{
	if (!isFileExist(filename))
		return true;

	QFile file(filename);
	file.setPermissions(QFile::WriteOwner);
	file.close();
	return file.remove();
}

bool BDFile::isFileExist(QString filepath)
{
	QFileInfo fi(filepath);
	return fi.isFile();
}

QList<QString> BDFile::getAllFiles(QString targetDir)
{
	QDir* dir = new QDir(targetDir);
	QStringList filter;
	filter = dir->entryList(QDir::Files);
	//		for (int i = 0; i<filter.size(); i++)
	//		{
	//			int set;
	//			set = filter[i].indexOf(".");
	//			filter[i] = filter[i].left(set);
	//		}
	return std::move(filter);
}

QList<QString> BDFile::getAllFilePaths(QString targetDir)
{
	if (!targetDir.endsWith("\\") || !targetDir.endsWith("/")) {
		targetDir += "\\";
	}

	QDir* dir = new QDir(targetDir);
	QStringList filter;
	filter = dir->entryList(QDir::Files);
	for (int i = 0; i < filter.size(); i++)
	{
		filter[i] = targetDir + filter[i];
	}
	return filter;
}

/**
* @brief 由文件名获取其所在的目录的路径
* @param filePath
* @return
*/
QString BDFile::getDirPath(QString filePath)
{
	QFileInfo fi(filePath);
	QString dir = fi.absoluteDir().path();
	if (!dir.endsWith("\\") || !dir.endsWith("/"))
		dir = dir + "/";
	return dir;
}

QString BDFile::getFileName(QString filePath)
{
	QFileInfo fi(filePath);
	return fi.baseName();
}

bool BDFile::waitForFile(QString fileName)
{
	while (!BDFile::isFileExist(fileName)) {
		//等待一段时间
		QTime t;
		t.start();
		while (t.elapsed() < 1000) {
			QCoreApplication::processEvents();
		}
	}
	return true;
}

bool BDFile::makeDir(QString dirName)
{
	QDir dir(dirName);
	if (!dir.exists()) {
		return dir.mkdir(dirName);
	}
	return true;
}

#else
bool BDFile::importData(const std::string filePath, BDMatrix& data)
{
	ifstream ifsReadIn;
	string strFile = filePath;//文件路径
	ifsReadIn.open(strFile);
	if (!ifsReadIn.is_open())
	{
		return false;
	}

	string buff;
	getline(ifsReadIn, buff);
	stringstream ss;
	ss << buff;
	int ncol = 0;
	while (!ss.eof())
	{
		string str;
		ss >> str;
		if (!str.empty())
			ncol++;
	}

	ifsReadIn.seekg(0, std::ios::beg);

	int iCountLine = 0;
	double dTemp = 0;
	while (!ifsReadIn.eof())
	{
		ifsReadIn >> dTemp;
		iCountLine++;
	}
	ifsReadIn.close();

	iCountLine /= ncol;
	data.Init(iCountLine, ncol);

	ifsReadIn.open(strFile);
	for (int i = 0; i < iCountLine; i++)
	{
		for (int j = 0; j < ncol; j++)
		{
			ifsReadIn >> data[i][j];
		}
	}
	ifsReadIn.close();
	return true;
}

#endif // inQT

/***********************************************************************
BDEnvironmentModel
***********************************************************************/

const double BDEarth::OMEAGA = 7.292115E-5;                              // 地球自转角速度(rad/s)
const double BDEarth::ro = 6371.0087714e3;                           // 地球平均半径(m)，用于圆球体模型
const double BDEarth::a = 6378.137e3;                               // 地球赤道半径(m)，认为赤道为圆形
const double BDEarth::alpha = (1 / 298.257222101);                        // 地球扁率
const double BDEarth::b = BDEarth::a * (1 - BDEarth::alpha);                // 经线椭圆短半径(m)
const double BDEarth::u = 3.986004418E14;                           // 地球引力常数u，万有引力系数和地球质量的乘积
const double BDEarth::J2 = 1.08263E-3;                               // 地球引力摄动J2项
const double BDEarth::J3 = -2.54e-6;                                 // 地球引力摄动J3项
const double BDEarth::J4 = -1.61e-6;                                 // 地球引力摄动J4项
const double BDEarth::J = 1.5 * BDEarth::J2;                           // 地球一阶扁率系数
const double BDEarth::e2 = 1 - BDEarth::b * BDEarth::b / (BDEarth::a * BDEarth::a);    // 地球偏心率常数的平方

const double BDEarth::ATMO_P0 = 101325.0;                                 // 标准大气压(pa)
const double BDEarth::ATMO_RHO0 = 1.225;                                    // 海平面大气密度(kg/m3)
const double BDEarth::ATMO_R0 = 6356.766;                                 // 大气常数
const double BDEarth::ATMO_AIRGAMMA = 1.4;                                      // 空气比热比
const double BDEarth::ATMO_RGAS = 287.14;                                   // 普适气体常数J/(kg*K)
const double BDEarth::ATMO_K = 273.15;                                   // 摄氏开氏温度互换
const double BDEarth::ATMO_STBO = 5.6697e-12;                               // 斯蒂芬-玻尔兹曼常数

BDEarth::BDEarth()
{}

BDEarth::~BDEarth()
{
}

double BDEarth::calGravity(double hM)
{
	double gn = 9.80665;
	double r = ATMO_R0 * 1000;
	double G_g = gn * pow(r / (r + hM), 2);
	return G_g;
}

double BDEarth::g0()
{
	return calGravity(0);
}

double BDEarth::calTK(double hM)
{
	if (hM < 0)
		hM = 0;

	double hKM = hM / 1000;
	double zKM = hKM / (1 + hKM / ATMO_R0); // 位势高度
	double W = 0;

	//标准大气模型
	if (0 <= hKM && hKM <= 11.0191)
	{
		W = 1 - zKM / 44.3308;
		return 288.15 * W;
	}
	else if (11.0191 < hKM && hKM <= 20.0631)
	{
		return 216.650;
	}
	else if (20.0631 < hKM && hKM <= 32.1619)
	{
		W = 1 + (zKM - 24.9021) / 221.552;
		return 221.552 * W;
	}
	else if (32.1619 < hKM && hKM <= 47.3501)
	{
		W = 1 + (zKM - 39.7499) / 89.4107;
		return 250.350 * W;
	}
	else if (47.3501 < hKM && hKM <= 51.4125)
	{
		W = exp((48.6252 - zKM) / 7.9223);
		return 270.650;
	}
	else if (51.4125 < hKM && hKM <= 71.8020)
	{
		W = 1 - (zKM - 59.4390) / 88.2218;
		return 247.021 * W;
	}
	else if (71.8020 < hKM && hKM <= 86.0000)
	{
		W = 1 - (zKM - 78.0303) / 100.2950;
		return 200.590 * W;
	}
	else if (86.0000 < hKM && hKM <= 91.0000)
	{
		W = exp((87.2848 - zKM) / 5.4700);
		return 186.8700;
	}
	return 0;
}

double BDEarth::calTC(double hM)
{
	return calTK(hM) - ATMO_K;
}

double BDEarth::calSonicVelocity(double hM)
{
	double TK = calTK(hM); // 开氏度
	double a = sqrt(ATMO_AIRGAMMA * ATMO_RGAS * TK);
	return a;
}

double BDEarth::calMa(double hM, double v)
{
	double a = calSonicVelocity(hM);
	if (a <= 0)
		return LARGENUM;
	return v / a;
}

double BDEarth::calDensity(double hM)
{
	if (hM < 0)
		hM = 0;

	double hKM = hM / 1000;
	double zKM = hKM / (1 + hKM / ATMO_R0); // 位势高度 km
	double W = 0;

	//标准大气模型
	if (0 <= hKM && hKM <= 11.0191)
	{
		W = 1 - zKM / 44.3308;
		return pow(W, 4.2559) * ATMO_RHO0;
	}
	else if (11.0191 < hKM && hKM <= 20.0631)
	{
		W = exp((14.9647 - zKM) / 6.3416);
		return 1.5898 * 0.1 * W * ATMO_RHO0;
	}
	else if (20.0631 < hKM && hKM <= 32.1619)
	{
		W = 1 + (zKM - 24.9021) / 221.552;
		return 3.2722 * 0.01 * pow(W, -35.1629) * ATMO_RHO0;
	}
	else if (32.1619 < hKM && hKM <= 47.3501)
	{
		W = 1 + (zKM - 39.7499) / 89.4107;
		return 3.2618 * 0.001 * pow(W, -13.2011) * ATMO_RHO0;
	}
	else if (47.3501 < hKM && hKM <= 51.4125)
	{
		W = exp((48.6252 - zKM) / 7.9223);
		return 9.4920 * 0.0001 * W * ATMO_RHO0;
	}
	else if (51.4125 < hKM && hKM <= 71.8020)
	{
		W = 1 - (zKM - 59.4390) / 88.2218;
		return 2.5280 * 0.0001 * pow(W, 11.2011) * ATMO_RHO0;
	}
	else if (71.8020 < hKM && hKM <= 86.0000)
	{
		W = 1 - (zKM - 78.0303) / 100.2950;
		return 1.7632 * 0.00001 * pow(W, 16.0816) * ATMO_RHO0;
	}
	else if (86.0000 < hKM && hKM <= 91.0000)
	{
		W = exp((87.2848 - zKM) / 5.4700);
		return 3.6411 * 0.000001 * W * ATMO_RHO0;
	}
	return 0;
}

double BDEarth::calAtmosHeight()
{
	return 91.0000e3;
}

double BDEarth::calStaticPressure(double hM)
{
	double rho = calDensity(hM);
	double T = calTK(hM);
	return 2.8705 *100 * rho * T;//modify*100
}

double BDEarth::calR0(double phi)
{
	double b = (1 - alpha);
	return a * b / sqrt(sin(phi) * sin(phi) + b * b * cos(phi) * cos(phi));
}

/***********************************************************************
BDOmega
***********************************************************************/

BDOmega::BDOmega(void)
{}
BDOmega::~BDOmega(void)
{}

void BDOmega::setLaunchInfo(const BDVector& info)
{
	launchInfo = info;
	double lambda0 = launchInfo[0];
	double B0 = launchInfo[1];
	double h0 = launchInfo[2];
	double A0 = launchInfo[3];

	auto X0 = BDCoordTransTool::LBH2Earth({ lambda0, B0, h0 });   // 发射点的地心系坐标
	Rx0 = BDCoordTransTool::Earth2Launchv(launchInfo, X0);      // 地心距矢量在发射系中的分量

	// 射向的正弦、余弦值
	double A0s = sin(A0);
	double A0c = cos(A0);

	//地球自转角速度在发射坐标系上的分量
	omeaga[0] = BDEarth::OMEAGA * cos(B0) * A0c;
	omeaga[1] = BDEarth::OMEAGA * sin(B0);
	omeaga[2] = -BDEarth::OMEAGA * cos(B0) * A0s;
	//牵连加速度在发射系下分量
	matrixA.Init(3, 3);
	matrixA[0][0] = pow(BDEarth::OMEAGA, 2.0) - pow(omeaga[0], 2.0);
	matrixA[1][1] = pow(BDEarth::OMEAGA, 2.0) - pow(omeaga[1], 2.0);
	matrixA[2][2] = pow(BDEarth::OMEAGA, 2.0) - pow(omeaga[2], 2.0);
	matrixA[0][1] = -omeaga[0] * omeaga[1];
	matrixA[1][0] = -omeaga[0] * omeaga[1];
	matrixA[0][2] = -omeaga[0] * omeaga[2];
	matrixA[2][0] = -omeaga[0] * omeaga[2];
	matrixA[1][2] = -omeaga[1] * omeaga[2];
	matrixA[2][1] = -omeaga[1] * omeaga[2];

	//柯式加速度在发射系下分量
	matrixB.Init(3, 3);
	matrixB[0][0] = 0;
	matrixB[0][1] = 2 * omeaga[2];
	matrixB[0][2] = -2 * omeaga[1];
	matrixB[1][0] = -2 * omeaga[2];
	matrixB[1][1] = 0;
	matrixB[1][2] = 2 * omeaga[0];
	matrixB[2][0] = 2 * omeaga[1];
	matrixB[2][1] = -2 * omeaga[0];
	matrixB[2][2] = 0;
}

BDPoint3 BDOmega::calGravity(const BDPoint3& p)
{
	BDPoint3 R = p + Rx0; // 发射系下的地心-导弹矢量
	double r = R.norm();  // 地心距

	auto LPH = BDCoordTransTool::Launch2LPH(launchInfo, p);
	double phi = LPH[1];  // 地心纬度

	// 引力加速度在地心矢径方向上的分量(ref.《远程火箭弹道学》Eq:2-2-17)
	double g_r = -BDEarth::u / pow(r, 2.0) * (1 + pow(BDEarth::a, 2.0) / pow(r, 2.0) * BDEarth::J * (1 - 5 * pow(sin(phi), 2.0)));

	// 引力加速度在地球自转轴上的分量(ref.《远程火箭弹道学》Eq:2-2-17)
	double g_omega = -2 * BDEarth::u * pow(BDEarth::a, 2.0) * BDEarth::J * sin(phi) / pow(r, 4.0);

	double gx = g_r * (p[0] + Rx0[0]) / r + g_omega * omeaga[0] / BDEarth::OMEAGA;
	double gy = g_r * (p[1] + Rx0[1]) / r + g_omega * omeaga[1] / BDEarth::OMEAGA;
	double gz = g_r * (p[2] + Rx0[2]) / r + g_omega * omeaga[2] / BDEarth::OMEAGA;
	return { gx, gy, gz };
}

BDPoint3 BDOmega::calAr(const BDPoint3& p)
{
	auto R = p + Rx0;      // 发射系下的地心-导弹矢量
	return matrixA * R;    // 发射系下的牵连加速度;
}

BDPoint3 BDOmega::calAc(const BDPoint3& v)
{
	return matrixB * v;    // 发射系下的柯式加速度
}

BDPoint3 BDOmega::calTotalGavity(const BDPoint3& p, const BDPoint3& v)
{
	auto G_g = calGravity(p);  // 发射系下的重力加速度
	auto ar = calAr(p);      // 发射系下的牵连加速度
	auto ac = calAc(v);      // 发射系下的柯式加速度
	return G_g + ar + ac;
}

double BDOmega::calThetaLocal(const BDPoint3& p, const BDPoint3& v, double theta)
{
	double vscalar = v.norm();
	if (iszero(vscalar))
		return theta;

	BDPoint3 R = p + Rx0;
	double r = R.norm();

	double s = (v * R) / r / vscalar;
	if (s < -1 + SMALLNUM)
		return -PI / 2;
	else if (s < 1 - SMALLNUM)
		return asin(s);
	else
		return PI / 2;
}

double BDOmega::getOmeagaz()
{
	return omeaga.z;
}

/***********************************************************************
BDInterpolation
***********************************************************************/

BDInterpolation::BDInterpolation()
{}
BDInterpolation::~BDInterpolation()
{
}
double BDInterpolation::interp1(int n_Xs, double* Xs, double* Data, double X)
{
	// Step1：确定Xs向位置
	int id;
	for (int i = 0; i != n_Xs - 1; i++)
	{
		id = i;
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		return Data[0] - ratio * (Data[1] - Data[0]);
	}
	else if (X >= Xs[n_Xs - 1])
	{
		ratio = (X - Xs[n_Xs - 1]) / (Xs[n_Xs - 1] - Xs[n_Xs - 2]);
		return Data[n_Xs - 1] + ratio * (Data[n_Xs - 1] - Data[n_Xs - 2]);
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		return Data[id] + (Data[id + 1] - Data[id]) * ratio;
	}
}

double BDInterpolation::interp1(int n_Xs, std::vector<double>& Xs, std::vector<double>& Data, double X)
{
	double tem;
	// Step1：确定Xs向位置
	int id = 0;
	for (int i = 0; i != n_Xs - 1; i++)
	{
		id = i;
		tem = X >= Xs[i];
		tem = Xs[i + 1];
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		return Data[0] - ratio * (Data[1] - Data[0]);
	}
	else if (X >= Xs[n_Xs - 1])
	{
		ratio = (X - Xs[n_Xs - 1]) / (Xs[n_Xs - 1] - Xs[n_Xs - 2]);
		return Data[n_Xs - 1] + ratio * (Data[n_Xs - 1] - Data[n_Xs - 2]);
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		return Data[id] + (Data[id + 1] - Data[id]) * ratio;
	}
}

double BDInterpolation::interp2(int n_Xs, double* Xs, int n_Ys, double* Ys, double* Data, double X, double Y)
{
	// Step1：确定Xs向位置
	int id = 0;
	for (int i = 0; i != n_Xs - 1; i++)
	{
		id = i;
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	double* Y_Data = new double[n_Ys];
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		for (int i = 0; i != n_Ys; i++)

			Y_Data[i] = Data[i * n_Xs] - ratio * (Data[i * n_Xs + 1] - Data[i * n_Xs]);
	}
	else if (X >= Xs[n_Xs - 1])
	{
		ratio = (X - Xs[n_Xs - 1]) / (Xs[n_Xs - 1] - Xs[n_Xs - 2]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + n_Xs - 1] + ratio * (Data[i * n_Xs + n_Xs - 1] - Data[i * n_Xs + n_Xs - 2]);
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + id] + (Data[i * n_Xs + id + 1] - Data[i * n_Xs + id]) * ratio;
	}

	// Step3：确定Ys向位置
	for (int i = 0; i != n_Ys; i++)
	{
		id = i;
		if (Y >= Ys[i] && Y <= Ys[i + 1])
			break;
		else
			continue;
	}

	// Step4：Ys向插值
	double	ResultVal;
	if (Y < Ys[0])
	{
		ratio = (Ys[0] - Y) / (Ys[1] - Ys[0]);
		ResultVal = Y_Data[0] - ratio * (Y_Data[1] - Y_Data[0]);
	}

	else if (Y >= Ys[n_Ys - 1])
	{
		ratio = (Y - Ys[n_Ys - 1]) / (Ys[n_Ys - 1] - Ys[n_Ys - 2]);
		ResultVal = Y_Data[n_Ys - 1] + ratio * (Y_Data[n_Ys - 1] - Y_Data[n_Ys - 2]);
	}
	else
	{
		ratio = (Y - Ys[id]) / (Ys[id + 1] - Ys[id]);
		ResultVal = Y_Data[id] + (Y_Data[id + 1] - Y_Data[id]) * ratio;
	}
	delete[] Y_Data;
	return ResultVal;
}

double BDInterpolation::interp3(int n_Xs, double* Xs, int n_Ys, double* Ys, int n_Zs, double* Zs, double* Data, double X, double Y, double Z)
{
	// Step1：Xs、Ys向二维插值
	double* Z_Data = new double[n_Zs];

	for (int i = 0; i != n_Zs; i++)
		Z_Data[i] = interp2(n_Xs, Xs, n_Ys, Ys, Data + i * n_Xs * n_Ys, X, Y);

	// Step2：确定Zs向位置
	int id;
	for (int i = 0; i != n_Zs; i++)
	{
		id = i;
		if (Z >= Zs[i] && Z <= Zs[i + 1])
			break;
		else
			continue;
	}

	// Step3：Zs向插值
	double ResultVal;
	double ratio;
	if (Z < Zs[0])
	{
		ratio = (Zs[0] - Z) / (Zs[1] - Zs[0]);
		ResultVal = Z_Data[0] - ratio * (Z_Data[1] - Z_Data[0]);
	}
	else if (Z >= Zs[n_Zs - 1])
	{
		ratio = (Z - Zs[n_Zs - 1]) / (Zs[n_Zs - 1] - Zs[n_Zs - 2]);
		ResultVal = Z_Data[n_Zs - 1] + ratio * (Z_Data[n_Zs - 1] - Z_Data[n_Zs - 2]);
	}

	else
	{
		ratio = (Z - Zs[id]) / (Zs[id + 1] - Zs[id]);
		ResultVal = Z_Data[id] + (Z_Data[id + 1] - Z_Data[id]) * ratio;
	}
	delete[] Z_Data;
	return ResultVal;
}

double BDInterpolation::interp4(int n_X1, double* X1s, int n_X2, double* X2s, int n_X3, double* X3s, int n_X4, double* X4s, double* Data, double X1, double X2, double X3, double X4)
{
	double w; //返回值

	double* Temp = new double[n_X4];
	int		Location; //三维差值表内数据的总个数

	for (int i = 0; i != n_X4; i++)
	{
		Location = n_X1 * n_X2 * n_X3 * i; //各组三维差值表起始位置相对Data的偏移量

		Temp[i] = interp3(n_X1, X1s, n_X2, X2s, n_X3, X3s, Data + Location, X1, X2, X3);//调用三维差值，计算给定第一、二、四维需差值数据，在第三维节点上的气动数据
	}
	w = interp1(n_X4, X4s, Temp, X4);
	delete[] Temp;
	return(w);
}

bool BDInterpolation::setData(BDMatrix m, int x_dim)
{
	std::vector<double> row;
	row.resize(m.GetNumColumns(), 0.0);

	std::vector< std::vector<double> > data;
	data.resize(m.GetNumRows(), row);

	for (int i = 0; i < m.GetNumRows(); i++)
		for (int j = 0; j < m.GetNumColumns(); j++)
			data[i][j] = m[i][j];
	return setData(data, x_dim);
}

bool BDInterpolation::setData(const std::vector< std::vector<double> >& interpolation_data, int x_dimen)
{
	this->interpolation_data = interpolation_data;
	this->x_dimen = x_dimen;
	this->total_dimen = interpolation_data[0].size();
	this->is_set_data = true;
	this->x_table = getXTable(x_dimen);

	//    mat = BDMatrix(interpolation_data);
	return true;
}

std::vector<std::vector<double> >& BDInterpolation::getData()
{
	return this->interpolation_data;
}

//BDMatrix BDInterpolation::getMatrix()
//{
//    return mat;
//}

std::vector<double> BDInterpolation::interp(double x)
{
	assert(x_dimen == 1);

	std::vector<double> x_data;
	x_data.push_back(x);

	return BDInterpolation::interp(x_data);
}

std::vector<double> BDInterpolation::interp(std::vector<double> x_data)
{
	assert(is_set_data);
	assert(x_dimen == x_data.size());

	//寻找附近的点
	std::vector<int> near_points = getNearPoints(x_data, x_table);

	//计算各维度的系数
	std::vector<double> k_list = getK(x_data, near_points, x_table);

	return getResult(k_list, near_points, x_table);
}

std::vector<int> BDInterpolation::getBinary(int num, int length)
{
	std::vector<int> res;
	res.resize(length);

	for (int i = length - 1; i >= 0; i--)
	{
		res[i] = num % 2;
		num /= 2;
	}
	return res;
}

std::vector<int> BDInterpolation::getNearPoints(std::vector<double>& target_x, std::vector<std::vector<double> >& x_table)
{
	std::vector<int> pos_list;
	int pos;
	double data;
	int size = target_x.size();
	pos_list.resize(size);

	for (int i = 0; i < size; i++)
	{
		pos = x_table[i].size() - 2;
		data = target_x[i];

		if (x_table[i].size() == 1)
		{
			pos = 0;
		}
		else
		{
			for (unsigned int j = 0; j < x_table[i].size() - 2; j++)
			{
				if (data < x_table[i][j + 1])
				{
					pos = j;
					break;
				}
			}
		}
		pos_list[i] = pos;
	}
	return pos_list;
}

std::vector<double> BDInterpolation::getK(std::vector<double>& target_x, std::vector<int>& component_pos, std::vector<std::vector<double> >& x_table)
{
	std::vector<double> k_list;
	double k;
	int length = target_x.size();
	k_list.resize(length);

	for (int i = 0; i < length; i++)
	{
		if (x_table[i].size() == 1)
		{
			k = 0;
			break;
		}

		int pos = component_pos[i];
		k = (target_x[i] - x_table[i][pos]) / (x_table[i][pos + 1] - x_table[i][pos]);
		k_list[i] = k;
	}

	return k_list;
}

std::vector<double> BDInterpolation::queryResultFromTable(std::vector<double>& x_data)
{
	std::vector<double> res;
	res.resize(total_dimen - x_dimen);

	int row_id = getRowId(x_data);
	if (row_id == -1)
		//找不到返回0
		return res;

	for (int i = 0; i < total_dimen - x_dimen; i++)
	{
		res[i] = interpolation_data[row_id][x_dimen + i];
	}
	return res;
}

std::vector<double> BDInterpolation::queryResultByComponentPosition(std::vector<unsigned int>& component_pos,
	std::vector<std::vector<double> >& x_table)
{
	std::vector<double> x_data;
	int length = component_pos.size();

	x_data.resize(length);

	for (int i = 0; i < length; i++)
	{
		if (component_pos[i] < x_table[i].size())
			x_data[i] = x_table[i][component_pos[i]];
		else
			x_data[i] = x_table[i][x_table[i].size() - 1];
	}

	return queryResultFromTable(x_data);
}

std::vector<double>  BDInterpolation::getResult(std::vector<double>& k_list,
	std::vector<int>& nearby_points,
	std::vector<std::vector<double> >& x_table)
{
	int dimen = k_list.size();
	int point_num = (int)(pow(2.0, dimen));
	int cols_num = total_dimen - x_dimen;

	std::vector<int> b_list;
	std::vector<double> res;
	std::vector<unsigned int> current_point_component_pos;

	res.resize(cols_num);

	current_point_component_pos.resize(nearby_points.size());

	for (int num = 0; num < point_num; num++)
	{
		b_list = getBinary(num, dimen);

		for (int i = 0; i < dimen; i++)
		{
			current_point_component_pos[i] = nearby_points[i] + b_list[i];
		}

		std::vector<double> y_list = queryResultByComponentPosition(current_point_component_pos, x_table);

		double k = 1;
		for (int i = 0; i < dimen; i++)
		{
			k *= 2 * b_list[i] * k_list[i] + 1 - b_list[i] - k_list[i];
		}

		for (int i = 0; i < cols_num; i++)
		{
			res[i] += k * y_list[i];
		}
	}

	return res;
}

//遍历总插值表
int BDInterpolation::getRowId(std::vector<double>& x_data)
{
	int row_id = -1;//没找到返回-1
	int total_row_num = interpolation_data.size();
	int x_data_dimen = x_data.size();
	bool is_this_row = true;

	for (int i = 0; i < total_row_num; i++)
	{
		is_this_row = true;

		for (int j = 0; j < x_data_dimen; j++)
		{
			if (interpolation_data[i][j] != x_data[j])
			{
				is_this_row = false;
				break;
			}
		}

		if (is_this_row)
		{
			row_id = i;
			return row_id;
		}
	}

	return row_id;
}

std::vector<std::vector<double> > BDInterpolation::getXTable(int x_dimen)
{
	std::vector< std::vector<double> > x_table;
	x_table.resize(x_dimen);

	int total_row_num = interpolation_data.size();

	for (int i = 0; i < x_dimen; i++)
	{
		for (int j = 0; j < total_row_num; j++)
		{
			if (!isContained(x_table[i], interpolation_data[j][i]))
			{
				int pos = getDataPosition(interpolation_data[j][i], x_table[i]);
				std::vector<double>::iterator iter = x_table[i].begin();
				iter += pos;
				x_table[i].insert(iter, 1, interpolation_data[j][i]);
			}
		}
	}
	return x_table;
}

/**
* 判断数组中是否包含某个数
*/
bool BDInterpolation::isContained(const std::vector<double>& list, double data)
{
	if (list.empty())
		return false;

	for (unsigned int i = 0; i < list.size(); i++)
	{
		if (list[i] == data)
			return true;
	}
	return false;
}

int BDInterpolation::getDataPosition(double data, std::vector<double>& list)
{
	int position = list.size();
	if (position == 0 || data >= list[position - 1])
		return position;

	for (unsigned int i = 0; i < list.size(); i++)
	{
		if (data < list[i])
		{
			position = i;
			break;
		}
	}
	return position;
}


/***********************************************************************
BDMath
***********************************************************************/
std::mt19937 BDMath::gen(time(nullptr));

BDMath::BDMath()
{}

BDMath::~BDMath()
{
}

void BDMath::bound(double minV, double& v, double maxV) {
	v = v > maxV ? maxV : v;
	v = v < minV ? minV : v;
}

#ifdef inQT
double BDMath::max(QVector<double> data) {
	if (data.size() == 0)
		return -1;

	double max = data[0];

	for (int i = 1; i < data.size(); i++) {
		if (max < data.at(i))
			max = data.at(i);
	}
	return max;
}

double BDMath::min(QVector<double> data) {
	if (data.size() == 0)
		return -1;

	double min = data[0];

	for (int i = 1; i < data.size(); i++) {
		if (min > data.at(i))
			min = data.at(i);
	}
	return min;
}

#endif

BDMatrix BDMath::uniform(double lower, double upper, int nrow, int ncol)
{
	std::uniform_real_distribution<> dis(lower, upper);
	BDMatrix m;
	m.Init(nrow, ncol);
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < ncol; j++)
		{
			m[i][j] = dis(gen);
		}
	}
	return m;
}

BDMatrix BDMath::normal(double mean, double std, int nrow, int ncol)
{
	std::normal_distribution<> dis(mean, std);
	BDMatrix m;
	m.Init(nrow, ncol);
	for (int i = 0; i < nrow; i++)
	{
		for (int j = 0; j < ncol; j++)
		{
			m[i][j] = dis(gen);
		}
	}
	return m;

	//        //        double rr = QRandomGenerator::global()->generate();
	//        double rr = rand();
	//        int m;
	//        double s = 65536.0;
	//        double ww = 2053.0;
	//        double vv = 13849.0;
	//        double tt = 0.0;
	//        for (int i = 1; i <= 12; i++)
	//        {
	//            rr = (rr)*ww + vv;
	//            m = (int)(rr / s);
	//            rr = rr - m*s;
	//            tt = tt + rr / s;
	//        }
	//        tt = mean + std*(tt - 6.0);
	//        if (tt > 6.0*std)
	//            tt = 6.0*std;
	//        return tt;
}

double BDMath::sign(double x)
{
	double z = 0.0;
	if (x > 0.0) z = 1.0;
	if (x < 0.0) z = -1.0;
	return z;
}


bool BDMath::solveEquation(std::function<double(double)> func, double lower, double upper, double& x)
{
	int maxloop = 20;           // 最大迭代次数
	double threshold = 1e-2;    // 阈值

	double para0 = lower;
	double para1 = upper;
	double parah = (para0 + para1) / 2.0;

	double y0 = func(para0);
	double y1 = func(para1);

	for (int kk = 0; kk < maxloop; kk++)
	{
		if (y0 * y1 > 0) {
			if (abs(y0) <= abs(y1))
				parah = para0;
			else
				parah = para1;
			break;
		}

		if (abs(y0) < threshold) {
			parah = para0;
			break;
		}

		if (abs(y1) < threshold) {
			parah = para1;
			break;
		}

		parah = (para0 + para1) / 2.0;

		double y = func(parah);

		if (y0 * y < 0) {
			para1 = parah;
			y1 = y;
		}
		else {
			para0 = parah;
			y0 = y;
		}
	}
	x = parah;
	if (abs(func(parah)) < threshold)
		return true;
	return false;
}

bool BDMath::solveEquation2(std::function<double(double)> func, double& x)
{
	int MaxIter = 7;
	double threshold = 1e-2;    // 阈值
	const int nOrder = 6;						// 连分式采用的最大阶数
	double IterX[nOrder], IterY[nOrder];		// 记录迭代过程值
	double IterB[nOrder];						// 记录迭代过程中的 b0、b1、b2……
	const int nMaxInitVal = 7;
	double SearchRange = 1000;					// 控制震荡搜索的空间

	// 外循环表示循环得到的结果不满足要求，需要以新的初值作为连分式的初值
	for (int IterOut = 0; IterOut < MaxIter; ++IterOut)
	{
		double InitX = x;		//内循环的初值
		if (IterOut > 0)
			InitX = IterX[nOrder - 1];

		int nInit = 0;
		double NoteInitVal = InitX;

		// 内循环表示连分式算法
		for (int IterIn = 0; IterIn < nOrder; ++IterIn)
		{
			/* 计算新的x值迭代 */
			if (IterIn < 2)
				IterX[IterIn] = InitX + 0.1 * IterIn;
			else
			{
				double TemNumera;		//临时分母项
				for (int i = IterIn - 1; i >= 0; --i)
				{
					if (i == IterIn - 1)
						TemNumera = IterB[i];
					else
					{
						TemNumera = IterB[i] - IterY[i] / TemNumera;
					}
				}
				IterX[IterIn] = TemNumera;
			}

			x = IterX[IterIn];

			IterY[IterIn] = func(IterX[IterIn]);

			/* 进行判断是否满足条件 */
			if (fabs(IterY[IterIn]) <= threshold)
			{
				// 满足条件返回;
				//x = IterX[IterIn];
				return true;
			}

			/* 根据新的x和y迭代出新的b */
			if (IterIn == 0)
				IterB[IterIn] = IterX[IterIn];
			else
			{
				// 迭代获得 bi 的过程
				double TemU = IterX[IterIn];
				for (int i = 0; i < IterIn; ++i)
				{
					// 非零检测
					if (fabs(TemU - IterB[i]) < 1e-8)
					{
						//自己写的，用于重新查找一个新的x值
						//TemU = IterB[i]+100
						TemU = 1e35;
					}
					else
						TemU = (IterY[IterIn] - IterY[i]) / (TemU - IterB[i]);
				}
				IterB[IterIn] = TemU;
			}

			if (IterIn != 0 && fabs(IterB[IterIn]) < 1e-8)
			{
				// 更换初值重新计算
				if (nInit < nMaxInitVal)
				{
					IterIn = -1;
					InitX = NoteInitVal + pow(double(-1), nInit) * ((nMaxInitVal - nInit) / 2 + 1) * 1000;
					nInit++;
				}
				else {
					// 随便给个值，是否收敛看运气
					IterB[IterIn] = 1;
				}
			}
		}
	}
	x = IterX[nOrder - 1];
	return false;
}

void BDMath::RungeKuta4(std::function<void(double, double*, double*)> func, double t, double step, double* y, unsigned int nDim)
{
	double cof[3];
	cof[0] = 0.5;
	cof[1] = 0.5;
	cof[2] = 1.0;

	double* y1 = new double[nDim];
	double* y2 = new double[nDim];
	double* dy = new double[nDim];

	func(t, y, dy);

	for (unsigned i = 0; i < nDim; i++)
		y2[i] = dy[i];

	for (unsigned k = 0; k < 3; k++)
	{
		for (unsigned i = 0; i < nDim; i++)
			y1[i] = y[i] + cof[k] * step * dy[i];

		func(t + cof[k] * step, y1, dy);

		for (unsigned i = 0; i < nDim; i++)
			y2[i] = y2[i] + dy[i] / cof[k];
	}

	for (unsigned i = 0; i < nDim; i++)
		y[i] = y[i] + y2[i] * step / 6.0;

	delete[] y1;
	delete[] y2;
	delete[] dy;
}

double BDMath::calAngle(double sinv, double cosv)
{
	if (cosv >= 0)
		return asin(sinv);
	else
	{
		if (sinv >= 0)
			return PI - asin(sinv);
		else
			return -PI - asin(sinv);
	}
}

double BDMath::calAngle2(double sinv, double cosv)
{
	if (cosv > 1.0 - SMALLNUM)
		return 0.0;
	else if (cosv < -1.0 + SMALLNUM)
		return PI;
	else
	{
		if (sinv >= 0.0)
			return acos(cosv);
		else
			return -acos(cosv);
	}
}


/***********************************************************************
BDTrajectoryTool
***********************************************************************/

BDTrajectoryTool::BDTrajectoryTool()
{
}

BDTrajectoryTool::~BDTrajectoryTool()
{
}

double BDTrajectoryTool::calT1(double N01)
{
	return sqrt(40.0 / (N01 - 1)); // 《远程火箭弹道学》Eq.9-2-3
}

double BDTrajectoryTool::calAlpha1(double t, double t1, double t2, double alphaMd)
{
	// 计算攻角规律的方法一，《远程火箭飞行动力学与制导》陈克俊，Eq.6-5-4
	if (t <= t1 || t > t2)
		return 0;

	double tm = t1 + (t2 - t1) / 5.0; // 三分之一处达到极值
	double k = (tm - t1) / (t2 - tm);
	double ft = PI * (t - t1) / (k * (t2 - t) + t - t1);
	double alphad = -fabs(alphaMd) * sin(ft) * sin(ft);
	return Deg2Rad(alphad);
}

double BDTrajectoryTool::calAlpha2(double t, double t1, double alphaMd, double a)
{
	// 计算攻角规律的方法二，《总体设计（上）》Eq.4.5-3
	double z = exp(a * (t1 - t));
	double alphad = -4 * fabs(alphaMd) * z * (z - 1);
	return Deg2Rad(alphad);
}

BDPoint3 BDTrajectoryTool::calGravity0(const BDPoint3& pos)
{
	double x = pos.x;
	double y = pos.y;
	double z = pos.z;
	double r = pos.norm();
	double gScalar = BDEarth::u / (r * r); // 引力加速度标量

	BDPoint3 g;
	g.x = -gScalar * (x / r);
	g.y = -gScalar * (y / r);
	g.z = -gScalar * (z / r);
	return g;
}

BDPoint3 BDTrajectoryTool::calGravity1(const BDPoint3& pos, const BDPoint3& vel)
{
	BDPoint3 g = calGravity0(pos);// 不考虑地球自转的引力

	double x = pos.x;
	double y = pos.y;
	double vx = vel.x;
	double vy = vel.y;
	BDPoint3 a1;
	a1.x = BDEarth::OMEAGA * BDEarth::OMEAGA * x + 2 * BDEarth::OMEAGA * vy;
	a1.y = BDEarth::OMEAGA * BDEarth::OMEAGA * y - 2 * BDEarth::OMEAGA * vx;
	a1.z = 0;
	return g + a1;
}

BDPoint3 BDTrajectoryTool::calGravity2(const BDPoint3& pos, const BDPoint3& vel)
{
	double x = pos.x;
	double y = pos.y;
	double z = pos.z;
	double vx = vel.x;
	double vy = vel.y;
	double vz = vel.z;
	double r = pos.norm();

	// 考虑J2项摄动的地球引力
	double r2 = pow(r, 2.0);
	double r3 = pow(r, 3.0);
	double J = 3 * BDEarth::J2 / 2;
	double fai = asin(pos.z / r);
	BDPoint3 g;
	g.x = -BDEarth::u * x / r3 * (1 + J * pow((BDEarth::a / r), 2.0) * (1 - 5 * pow(sin(fai), 2.0)));
	g.y = -BDEarth::u * y / r3 * (1 + J * pow((BDEarth::a / r), 2.0) * (1 - 5 * pow(sin(fai), 2.0)));
	g.z = -BDEarth::u * z / r3 * (1 + J * pow((BDEarth::a / r), 2.0) * (1 - 5 * pow(sin(fai), 2.0))) - 2 * BDEarth::u / r2 * J * pow(BDEarth::a / r, 2.0) * sin(fai);

	// 计算地固系总加速度
	double w = BDEarth::OMEAGA * BDEarth::OMEAGA;
	BDPoint3 acc;
	acc.x = g.x + w * x + 2 * BDEarth::OMEAGA * vx;
	acc.y = g.y + w * y - 2 * BDEarth::OMEAGA * vy;
	acc.z = g.z;
	return acc;
}

//    void BDTrajectoryTool::freeTraj(
//            const BDVector &launchInfo,
//            const BDVector &initTPVM,
//            double tspan,
//            BDMatrix *resultTP)
//    {
//        BDTrajSolver s;
//        s.setInterfaces(launchInfo);
//        s.setInitData(initTPVM);
//        s.setExportFunction([](BDTrajSolver::BDTrajectoryStateParam &param, BDVector &stepdata)
//        {
//            stepdata[0] = param.t;
//            stepdata[1] = param.X[0];
//            stepdata[2] = param.X[1];
//            stepdata[3] = param.X[2];
//        },
//        [](std::vector<std::string> &title)
//        {
//            title.resize(4);
//            title[0] = "t/s";
//            title[1] = "x/m";
//            title[2] = "y/m";
//            title[3] = "z/m";
//        });

//        double step = 1.0; // 仿真计算时间步长
//        double eps = TIME_STEP/100.0;

//        bool isContinue = true;
//        double currentTime = initTPVM[0];
//        double stop = currentTime + tspan;
//        for(; currentTime + eps + step < stop && isContinue; currentTime += step)
//            isContinue = s.oneStep(step);

//        s.oneStep(stop - currentTime);
//        s.getResult(resultTP);
//    }

BDMatrix BDTrajectoryTool::predictBallisticTraj(const BDPoint3& pos, const BDPoint3& vel, double endT)
{
	if (endT < 0)
		endT = LARGENUM;

	double t_predict = 0.0;
	auto p_predict = pos;
	auto v_predict = vel;
	auto a_predict = calGravity2(p_predict, v_predict);

	std::vector<std::vector<double> > result;
	result.push_back({ t_predict,
					  p_predict.x, p_predict.y, p_predict.z,
					  v_predict.x, v_predict.y, v_predict.z,
					  a_predict.x, a_predict.y, a_predict.z });

	double deltaT = 1;
	double H_predict = LARGENUM;
	while (t_predict < endT)
		//while (t_predict < endT && H_predict > 80e3)
	{
		t_predict = t_predict + deltaT;
		p_predict = p_predict + v_predict * deltaT;
		v_predict = v_predict + a_predict * deltaT;
		a_predict = calGravity2(p_predict, v_predict);

		auto LLA_predict = BDCoordTransTool::Earth2LPH(p_predict);
		H_predict = LLA_predict[2];

		result.push_back({ t_predict,
						  p_predict.x, p_predict.y, p_predict.z,
						  v_predict.x, v_predict.y, v_predict.z,
						  a_predict.x, a_predict.y, a_predict.z });
	}

	return BDMatrix(result);
}


/***********************************************************************
BDTrajectoryTool
***********************************************************************/

BDVector::BDVector()
{}

BDVector::BDVector(const std::vector<double>& list)
{
	resize(list.size());
	for (unsigned i = 0; i < list.size(); i++) {
		(*this)[i] = list.at(i);
	}
}

BDVector::~BDVector()
{}

void BDVector::setElement(unsigned int idx, double value)
{
	assert(idx < this->getLength());
	(*this)[idx] = value;
}

double BDVector::getElement(unsigned int idx) const
{
	assert(idx < this->getLength());
	return (*this)[idx];
}

double BDVector::firstElement()
{
	return getElement(0);
}

double BDVector::lastElement()
{
	return getElement(getLength() - 1);
}

unsigned int BDVector::getLength() const
{
	return size();
}

#ifdef inQT
void BDVector::print() const
{
	auto ncol = getLength();
	QString str;
	for (unsigned i = 0; i < ncol - 1; i++)
	{
		str += QString::number((*this)[i]) + QString("  ");
	}
	str += QString::number((*this)[ncol - 1]);
	BDLOG str;
}
#endif

double BDVector::max() const
{
	if (size() == 0)
		return 0;
	double max = this->at(0);
	for (int i = 1; i < size(); i++)
	{
		if (max < this->at(i))
			max = this->at(i);
	}
	return max;
}

void BDVector::forEachElement(std::function<void(unsigned int, double&)> op)
{
	for (unsigned int idx = 0; idx < getLength(); idx++)
	{
		op(idx, (*this)[idx]);
	}
}

BDEulerAngle::BDEulerAngle() {
}

BDEulerAngle::~BDEulerAngle()
{
}

double BDEulerAngle::calEta(double alpha, double beta)
{
	// 《远程火箭弹道学》, Eq.(6-1-11)
	return acos(cos(alpha) * cos(beta));
}

BDPoint3 BDEulerAngle::calThetaSigmaGammav(const BDPoint3& v)
{
	// 《弹道导弹弹道学》, P25
	double vx = v.x, vy = v.y, vz = v.z;

	//计算发射系弹道倾角
	double theta;
	if (vx < -SMALLNUM)
	{
		if (vy > 0)
			theta = PI - atan(fabs(vy / vx));
		else
			theta = -PI + atan(fabs(vy / vx));
	}
	else if (vx < SMALLNUM)
	{
		if (vy >= 0)
			theta = PI / 2 - SMALLNUM;//垂直发射，初始值
		else
			theta = -PI / 2 + SMALLNUM;
	}
	else
	{
		theta = atan(vy / vx);
	}

	//计算发射系弹道偏角
	double vel = v.norm();
	double sigma;
	if (vel < SMALLNUM)
		sigma = 0;
	else
		sigma = -asin(vz / vel); // 弹道偏角，以左为正

	double gammav = 0;
	return { theta, sigma, gammav };
}

BDPoint3 BDEulerAngle::calPhiPsiGamma(BDPoint3 theta_sigma_gammav, BDPoint3 alpha_beta)
{
	double theta = theta_sigma_gammav.x;
	double sigma = theta_sigma_gammav.y;
	double gammav = theta_sigma_gammav.z;
	double alpha = alpha_beta.x;
	double beta = alpha_beta.y;

	double psi = asin(cos(alpha) * cos(beta) * sin(sigma) - sin(alpha) * cos(sigma) * sin(gammav) + cos(alpha) * sin(beta) * cos(sigma) * sin(gammav));

	double fais = (cos(alpha) * cos(beta) * sin(theta) * cos(sigma) + sin(alpha) * sin(theta) * sin(sigma) * sin(gammav) + sin(alpha) * cos(theta) * cos(gammav) - cos(alpha) * sin(beta) * sin(theta) * sin(sigma) * cos(gammav) + cos(alpha) * sin(beta) * cos(theta) * sin(gammav)) / cos(psi);
	double faic = (cos(alpha) * cos(beta) * cos(theta) * cos(sigma) + sin(alpha) * cos(theta) * sin(sigma) * sin(gammav) - sin(alpha) * sin(theta) * cos(gammav) - cos(alpha) * sin(beta) * cos(theta) * sin(sigma) * cos(gammav) - cos(alpha) * sin(beta) * sin(theta) * sin(gammav)) / cos(psi);
	if ((faic * faic + fais * fais - 1 >= 0.1) || ((faic * faic + fais * fais - 1 <= -0.1)))
	{
		cout << "角度计算错误" << endl;
	}

	double fai;
	if (fais > 1)
		fais = 1 - SMALLNUM;
	if (fais < -1)
		fais = -1 + SMALLNUM;

	if (faic >= 0)
		fai = asin(fais);
	else
	{
		if (fais >= 0)
			fai = PI - asin(fais);
		else
			fai = -PI - asin(fais);
	}
	return { fai, psi, 0 };
}

BDPoint3 BDEulerAngle::calPhiPsiGamma(double B0, double A0, double t, const BDPoint3& phi_psi_gamma_a)
{
	double X = cos(B0) * cos(A0);
	double Y = sin(B0);
	double Z = -cos(B0) * sin(A0);
	double S = sin(BDEarth::OMEAGA * t);
	double C = cos(BDEarth::OMEAGA * t);
	BDMatrix m1(3, 3);           // 发射到发惯系
	m1[0][0] = X * X * (1 - C) + C;
	m1[0][1] = X * Y * (1 - C) - Z * S;
	m1[0][2] = X * Z * (1 - C) + Y * S;
	m1[1][0] = X * Y * (1 - C) + Z * S;
	m1[1][1] = Y * Y * (1 - C) + C;
	m1[1][2] = Y * Z * (1 - C) - X * S;
	m1[2][0] = X * Z * (1 - C) - Y * S;
	m1[2][1] = Y * Z * (1 - C) + X * S;
	m1[2][2] = Z * Z * (1 - C) + C;

	// 发惯系姿态角
	double phia = phi_psi_gamma_a[0];
	double psia = phi_psi_gamma_a[1];
	double gammaa = phi_psi_gamma_a[2];
	BDMatrix m2(3, 3);                 // 发惯到弹体系转换矩阵
	m2[0][0] = cos(phia) * cos(psia);
	m2[0][1] = sin(phia) * cos(psia);
	m2[0][2] = -sin(psia);
	m2[1][0] = cos(phia) * sin(psia) * sin(gammaa) - sin(phia) * cos(gammaa);
	m2[1][1] = sin(phia) * sin(psia) * sin(gammaa) + cos(phia) * cos(gammaa);
	m2[1][2] = cos(psia) * sin(gammaa);
	m2[2][0] = cos(phia) * sin(psia) * cos(gammaa) + sin(phia) * sin(gammaa);
	m2[2][1] = sin(phia) * sin(psia) * cos(gammaa) - cos(phia) * sin(gammaa);
	m2[2][2] = cos(psia) * cos(gammaa);

	auto m = m2 * m1; // 发射到弹体系转换矩阵

	// 计算发射系姿态角
	double psi = asin(-m[0][2]);
	double sin_phi = m[0][1] / cos(psi);
	double cos_phi = m[0][0] / cos(psi);
	double sin_gamma = m[1][2] / cos(psi);
	double cos_gamma = m[2][2] / cos(psi);
	double phi = BDMath::calAngle2(sin_phi, cos_phi);
	double gamma = BDMath::calAngle2(sin_gamma, cos_gamma);
	return BDPoint3(phi, psi, gamma);
}

BDPoint3 BDEulerAngle::calPhiPsiGammaa(double B0, double A0, double t, const BDPoint3& phi_psi_gamma)
{
	double X = cos(B0) * cos(A0);
	double Y = sin(B0);
	double Z = -cos(B0) * sin(A0);
	double S = sin(BDEarth::OMEAGA * t);
	double C = cos(BDEarth::OMEAGA * t);
	BDMatrix m1(3, 3);
	m1[0][0] = X * X * (1 - C) + C;
	m1[0][1] = X * Y * (1 - C) - Z * S;
	m1[0][2] = X * Z * (1 - C) + Y * S;
	m1[1][0] = X * Y * (1 - C) + Z * S;
	m1[1][1] = Y * Y * (1 - C) + C;
	m1[1][2] = Y * Z * (1 - C) - X * S;
	m1[2][0] = X * Z * (1 - C) - Y * S;
	m1[2][1] = Y * Z * (1 - C) + X * S;
	m1[2][2] = Z * Z * (1 - C) + C;
	m1 = m1.Transpose(); // 发惯系到发射系

	double phi = phi_psi_gamma[0];
	double psi = phi_psi_gamma[1];
	double gamma = phi_psi_gamma[2];
	BDMatrix m2(3, 3);

	//        m2.SetElement(0, 0, cos(phi)*cos(psi));
	//        m2.SetElement(0, 1, cos(phi)*sin(psi)*sin(gamma) - sin(phi)*cos(gamma));
	//        m2.SetElement(0, 2, cos(phi)*sin(psi)*cos(gamma) + sin(phi)*sin(gamma));
	//        m2.SetElement(1, 0, sin(phi)*cos(psi));
	//        m2.SetElement(1, 1, sin(phi)*sin(psi)*sin(gamma) + cos(phi)*cos(gamma));
	//        m2.SetElement(1, 2, sin(phi)*sin(psi)*cos(gamma) - cos(phi)*sin(gamma));
	//        m2.SetElement(2, 0, -sin(psi));
	//        m2.SetElement(2, 1, cos(psi)*sin(gamma));
	//        m2.SetElement(2, 2, cos(psi)*cos(gamma));

	double c_phi = cos(phi);
	double s_phi = sin(phi);
	double c_psi = cos(psi);
	double s_psi = sin(psi);
	double s_gamma = sin(gamma);
	double c_gamma = cos(gamma);
	m2.SetElement(0, 0, c_phi * c_psi);
	m2.SetElement(0, 1, c_phi * s_psi * s_gamma - s_phi * c_gamma);
	m2.SetElement(0, 2, c_phi * s_psi * c_gamma + s_phi * s_gamma);
	m2.SetElement(1, 0, s_phi * c_psi);
	m2.SetElement(1, 1, s_phi * s_psi * s_gamma + c_phi * c_gamma);
	m2.SetElement(1, 2, s_phi * s_psi * c_gamma - c_phi * s_gamma);
	m2.SetElement(2, 0, -s_psi);
	m2.SetElement(2, 1, c_psi * s_gamma);
	m2.SetElement(2, 2, c_psi * c_gamma);

	m2 = m2.Transpose(); // 发射系到弹体系
	auto m = m2 * m1; // 发惯系到弹体系

	// 计算发惯系姿态角
	double psia = asin(-m[0][2]);
	double sin_phia = m[0][1] / cos(psia);
	double cos_phia = m[0][0] / cos(psia);
	double sin_gammaa = m[1][2] / cos(psia);
	double cos_gammaa = m[2][2] / cos(psia);
	double phia = BDMath::calAngle2(sin_phia, cos_phia);
	double gammaa = BDMath::calAngle2(sin_gammaa, cos_gammaa);
	return BDPoint3(phia, psia, gammaa);
}

BDPoint3 BDEulerAngle::calAlphaBetaNu(BDPoint3 theta_sigma_gammav, BDPoint3 phi_psi_gamma)
{
	// 《远程火箭飞行动力学与制导》陈克俊.Eq.3-2-32
	double theta = theta_sigma_gammav[0];
	double sigma = theta_sigma_gammav[1];
	double phi = phi_psi_gamma[0];
	double psi = phi_psi_gamma[1];
	double gamma = phi_psi_gamma[2];

	double sin_beta = cos(phi - theta) * cos(sigma) * sin(psi) * cos(gamma)
		+ sin(phi - theta) * cos(sigma) * sin(gamma)
		- sin(sigma) * cos(psi) * cos(gamma);
	double beta = asin(sin_beta);

	double cos_alpha = (cos(sigma) * cos(theta - phi) * cos(psi) + sin(sigma) * sin(psi)) / cos(beta);
	double sin_alpha = -(-cos(psi) * sin(gamma) * sin(sigma) + cos(sigma) * (cos(gamma) * sin(theta - phi) + cos(theta - phi) * sin(gamma) * sin(psi))) / cos(beta);
	double alpha = BDMath::calAngle2(sin_alpha, cos_alpha);

	double sin_nu = (cos(alpha) * cos(psi) * sin(gamma) - sin(alpha) * sin(psi)) / cos(sigma);
	double cos_nu = (cos(beta) * cos(psi) * cos(gamma) + sin(beta) * (sin(gamma) * cos(psi) * sin(alpha) + cos(alpha) * sin(psi))) / cos(sigma);
	double nu = BDMath::calAngle2(sin_nu, cos_nu);

	return { alpha, beta, nu };
}



/******************************************************************************************************/
//C_Matrix33
/******************************************************************************************************/

C_Matrix33::C_Matrix33()
{
	a11 = 0;
	a12 = 0;
	a13 = 0;
	a21 = 0;
	a22 = 0;
	a23 = 0;
	a31 = 0;
	a32 = 0;
	a33 = 0;
}

C_Matrix33::C_Matrix33(double a_set11, double a_set12, double a_set13, double a_set21, double a_set22, double a_set23, double a_set31, double a_set32, double a_set33)
{
	a11 = a_set11;
	a12 = a_set12;
	a13 = a_set13;
	a21 = a_set21;
	a22 = a_set22;
	a23 = a_set23;
	a31 = a_set31;
	a32 = a_set32;
	a33 = a_set33;
}
C_Matrix33::~C_Matrix33()
{
	a11 = 0;
	a12 = 0;
	a13 = 0;
	a21 = 0;
	a22 = 0;
	a23 = 0;
	a31 = 0;
	a32 = 0;
	a33 = 0;
}

void C_Matrix33::fn_SetElems(double a_set11, double a_set12, double a_set13, double a_set21, double a_set22, double a_set23, double a_set31, double a_set32, double a_set33)
{
	a11 = a_set11;
	a12 = a_set12;
	a13 = a_set13;
	a21 = a_set21;
	a22 = a_set22;
	a23 = a_set23;
	a31 = a_set31;
	a32 = a_set32;
	a33 = a_set33;
}

double C_Matrix33::fn_CalDet()
{
	double det = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31 - a12 * a21 * a33 - a11 * a23 * a32;
	return det;
}

C_Matrix33 C_Matrix33::fn_CalInv(C_Matrix33 IN_CM33)
{
	C_Matrix33 inv;
	double det = IN_CM33.a11 * IN_CM33.a22 * IN_CM33.a33 + IN_CM33.a12 * IN_CM33.a23 * IN_CM33.a31 + IN_CM33.a13 * IN_CM33.a21 * IN_CM33.a32
		- IN_CM33.a13 * IN_CM33.a22 * IN_CM33.a31 - IN_CM33.a12 * IN_CM33.a21 * IN_CM33.a33 - IN_CM33.a11 * IN_CM33.a23 * IN_CM33.a32;

	if (det == 0)
	{
		cout << "ERROR" << endl;
		system("pause");
		exit(0);
	}

	inv.a11 = (IN_CM33.a22 * IN_CM33.a33 - IN_CM33.a32 * IN_CM33.a23) / det;
	inv.a12 = (IN_CM33.a32 * IN_CM33.a13 - IN_CM33.a12 * IN_CM33.a33) / det;
	inv.a13 = (IN_CM33.a12 * IN_CM33.a23 - IN_CM33.a13 * IN_CM33.a22) / det;
	inv.a11 = (IN_CM33.a22 * IN_CM33.a33 - IN_CM33.a32 * IN_CM33.a23) / det;
	inv.a12 = (IN_CM33.a32 * IN_CM33.a13 - IN_CM33.a12 * IN_CM33.a33) / det;
	inv.a13 = (IN_CM33.a12 * IN_CM33.a23 - IN_CM33.a13 * IN_CM33.a22) / det;
	inv.a21 = (IN_CM33.a23 * IN_CM33.a31 - IN_CM33.a33 * IN_CM33.a21) / det;
	inv.a22 = (IN_CM33.a33 * IN_CM33.a11 - IN_CM33.a13 * IN_CM33.a31) / det;
	inv.a23 = (IN_CM33.a13 * IN_CM33.a21 - IN_CM33.a11 * IN_CM33.a23) / det;
	inv.a31 = (IN_CM33.a21 * IN_CM33.a32 - IN_CM33.a31 * IN_CM33.a22) / det;
	inv.a32 = (IN_CM33.a31 * IN_CM33.a12 - IN_CM33.a11 * IN_CM33.a32) / det;
	inv.a33 = (IN_CM33.a11 * IN_CM33.a22 - IN_CM33.a12 * IN_CM33.a21) / det;
	return inv;
};

	//求转置函数
C_Matrix33 C_Matrix33::fn_CalTran(C_Matrix33 IN_CM33)
{
	C_Matrix33 Tran;
	Tran.a11 = IN_CM33.a11;
	Tran.a12 = IN_CM33.a21;
	Tran.a13 = IN_CM33.a31;
	Tran.a21 = IN_CM33.a12;
	Tran.a22 = IN_CM33.a22;
	Tran.a23 = IN_CM33.a32;
	Tran.a31 = IN_CM33.a13;
	Tran.a32 = IN_CM33.a23;
	Tran.a33 = IN_CM33.a33;

	return Tran;
};

C_Matrix33 operator*(C_Matrix33 IN_CM33_L, C_Matrix33 IN_CM33_R)
{
	C_Matrix33 result;
	result.a11 = IN_CM33_L.a11 * IN_CM33_R.a11 + IN_CM33_L.a12 * IN_CM33_R.a21 + IN_CM33_L.a13 * IN_CM33_R.a31;
	result.a12 = IN_CM33_L.a11 * IN_CM33_R.a12 + IN_CM33_L.a12 * IN_CM33_R.a22 + IN_CM33_L.a13 * IN_CM33_R.a32;
	result.a13 = IN_CM33_L.a11 * IN_CM33_R.a13 + IN_CM33_L.a12 * IN_CM33_R.a23 + IN_CM33_L.a13 * IN_CM33_R.a33;
	result.a21 = IN_CM33_L.a21 * IN_CM33_R.a11 + IN_CM33_L.a22 * IN_CM33_R.a21 + IN_CM33_L.a23 * IN_CM33_R.a31;
	result.a22 = IN_CM33_L.a21 * IN_CM33_R.a12 + IN_CM33_L.a22 * IN_CM33_R.a22 + IN_CM33_L.a23 * IN_CM33_R.a32;
	result.a23 = IN_CM33_L.a21 * IN_CM33_R.a13 + IN_CM33_L.a22 * IN_CM33_R.a23 + IN_CM33_L.a23 * IN_CM33_R.a33;
	result.a31 = IN_CM33_L.a31 * IN_CM33_R.a11 + IN_CM33_L.a32 * IN_CM33_R.a21 + IN_CM33_L.a33 * IN_CM33_R.a31;
	result.a22 = IN_CM33_L.a31 * IN_CM33_R.a12 + IN_CM33_L.a32 * IN_CM33_R.a22 + IN_CM33_L.a33 * IN_CM33_R.a32;
	result.a23 = IN_CM33_L.a31 * IN_CM33_R.a13 + IN_CM33_L.a32 * IN_CM33_R.a23 + IN_CM33_L.a33 * IN_CM33_R.a33;
	return result;

};

void C_Matrix33::fn_Disp()
{
	cout << "[" << a11 << "," << a12 << "," << a13 << "\n"
		<< " " << a21 << "," << a22 << "," << a23 << "\n"
		<< " " << a31 << "," << a32 << "," << a33 << "]" << endl;
};


//***************************************************************************************
C_Vector3::C_Vector3()
{
	x = 0;
	y = 0;
	z = 0;
};

C_Vector3::C_Vector3(double x_in, double y_in, double z_in)
{
	x = x_in;
	y = y_in;
	z = z_in;
}
C_Vector3::C_Vector3(const BDPoint3& in_arr)
{
	x = in_arr.x;
	y = in_arr.y;
	z = in_arr.z;
}
C_Vector3::C_Vector3(const C_Vector3& in_arr)
{
	x = in_arr.x;
	y = in_arr.y;
	z = in_arr.z;
}
C_Vector3::C_Vector3(const double in_arr[3])
{
	x = in_arr[0];
	y = in_arr[1];
	z = in_arr[2];
}
C_Vector3::C_Vector3(vector<double> in_arr)
{
	if (in_arr.size() != 3)
	{
		assert(1);
	}
	else
	{
		x = in_arr[0];
		y = in_arr[1];
		z = in_arr[2];
	}
}


C_Vector3::~C_Vector3()
{
}

//设置向量值函数
void C_Vector3::SetElems(double Xset, double Yset, double Zset)
{
	x = Xset;
	y = Yset;
	z = Zset;
}

//向量求和函数
C_Vector3 operator+(const C_Vector3& CV3_1, const C_Vector3& CV3_2)
{
	C_Vector3 result;
	result.x = CV3_1.x + CV3_2.x;
	result.y = CV3_1.y + CV3_2.y;
	result.z = CV3_1.z + CV3_2.z;
	return result;
}

//向量求差函数
C_Vector3 operator-(const C_Vector3& CV3_1, const C_Vector3& CV3_2)
{
	C_Vector3 result;
	result.x = CV3_1.x - CV3_2.x;
	result.y = CV3_1.y - CV3_2.y;
	result.z = CV3_1.z - CV3_2.z;
	return result;
}

//矩阵*向量函数
C_Vector3 operator*(const C_Matrix33& CV33, const C_Vector3& CV3)
{
	C_Vector3 result;
	result.x = CV33.a11 * CV3.x + CV33.a12 * CV3.y + CV33.a13 * CV3.z;
	result.y = CV33.a21 * CV3.x + CV33.a22 * CV3.y + CV33.a23 * CV3.z;
	result.z = CV33.a31 * CV3.x + CV33.a32 * CV3.y + CV33.a33 * CV3.z;
	return result;
}

//窗口显示向量
void C_Vector3::Disp()
{
	cout << fixed << setprecision(3) << "(" << x << ", " << y << ", " << z << ")" << endl;
};

//向量*一个数，返回值为一个向量
C_Vector3 operator*(const double& IN_num, const C_Vector3& IN_CV3)
{
	C_Vector3 result;
	result.x = IN_num * IN_CV3.x;
	result.y = IN_num * IN_CV3.y;
	result.z = IN_num * IN_CV3.z;
	return result;
}

//向量*一个数，返回值为一个向量
C_Vector3 operator*(const C_Vector3& IN_CV3, const double& IN_num)
{
	C_Vector3 result;
	result.x = IN_num * IN_CV3.x;
	result.y = IN_num * IN_CV3.y;
	result.z = IN_num * IN_CV3.z;
	return result;
}

C_Vector3 operator/(const C_Vector3& IN_CV3, const double& IN_num)
{
	C_Vector3 result;
	result.x = IN_CV3.x / IN_num;
	result.y = IN_CV3.y / IN_num;
	result.z = IN_CV3.z / IN_num;

	return result;
}

C_Vector3 operator*(C_Vector3 IN_CV3_L, C_Vector3 IN_CV3_R)
{
	C_Vector3 result;
	result.x = IN_CV3_L.y * IN_CV3_R.z - IN_CV3_L.z * IN_CV3_R.y;
	result.y = IN_CV3_L.z * IN_CV3_R.x - IN_CV3_L.x * IN_CV3_R.z;
	result.z = IN_CV3_L.x * IN_CV3_R.y - IN_CV3_L.y * IN_CV3_R.x;
	return result;
}

double C_Vector3::Mod() const
{
	double ModLength = sqrt(x * x + y * y + z * z);
	return ModLength;
}

double C_Vector3::Distance(C_Vector3 in_arr)
{
	return sqrt(pow((x - in_arr.x), 2) + pow((y - in_arr.y), 2) + pow((z - in_arr.z), 2));
}

double C_Vector3::innerProduct(const C_Vector3& other)
{
	return x * other.x + y * other.y + z * other.z;
}

C_Vector3 C_Vector3::Direction() const
{
	double ModLength = sqrt(x * x + y * y + z * z);

	C_Vector3 DirectionVector;
	DirectionVector.x = x / ModLength;
	DirectionVector.y = y / ModLength;
	DirectionVector.z = z / ModLength;

	return DirectionVector;
}

C_Vector3 C_Vector3::operator+=(const C_Vector3& IN_CV3)
{
	(*this) = (*this) + IN_CV3;
	return (*this);
}

C_Vector3 C_Vector3::operator-=(const C_Vector3& IN_CV3)
{
	(*this) = (*this) - IN_CV3;
	return (*this);
}

C_Vector3 C_Vector3::operator*=(const double IN_num)
{
	(*this) = (*this) * IN_num;
	return (*this);
}

C_Vector3 C_Vector3::operator*=(const C_Matrix33& IN_CV33)
{
	(*this) = IN_CV33 * (*this);
	return (*this);
}

C_Vector3 C_Vector3::operator/=(const double IN_num)
{
	(*this) = (*this) / IN_num;
	return (*this);
}

double C_Vector3::operator[](int index)
{
	switch (index)
	{
	case 0:
		return x;
	case 1:
		return y;
	case 2:
		return z;
	default:
		assert(false);
		return x;
	}
}

bool C_Vector3::operator==(C_Vector3 point)
{
	if (x == point.x && y == point.y && z == point.z)
		return true;
	else
		return false;
}

bool C_Vector3::operator!=(C_Vector3 point)
{
	if (x == point.x && y == point.y && z == point.z)
		return false;
	else
		return true;
}

C_Vector3 C_Vector3::operator-()
{
	return C_Vector3(-x,-y,-z);
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//************函数库部分*************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////坐标转换函数/////////////////////////////////////////////////////////

//发射系到弹体系，航天一院所提供
C_Matrix33 getMatrix_LaunchToBody(double IN_dPhi, double IN_dPsi, double IN_dGamma)
{
	C_Matrix33 LaunchToBody;
	double cPhi = cos(IN_dPhi);
	double sPhi = sin(IN_dPhi);
	double cPsi = cos(IN_dPsi);
	double sPsi = sin(IN_dPsi);
	double cGamma = cos(IN_dGamma);
	double sGamma = sin(IN_dGamma);
	LaunchToBody.a11 = cPhi * cPsi;
	LaunchToBody.a12 = sPhi * cPsi;
	LaunchToBody.a13 = -sPsi;
	LaunchToBody.a21 = cPhi * sPsi * sGamma - sPhi * cGamma;
	LaunchToBody.a22 = sPhi * sPsi * sGamma + cPhi * cGamma;
	LaunchToBody.a23 = cPsi * sGamma;
	LaunchToBody.a31 = cPhi * sPsi * cGamma + sPhi * sGamma;
	LaunchToBody.a32 = sPhi * sPsi * cGamma - cPhi * sGamma;
	LaunchToBody.a33 = cPsi * cGamma;
	//LaunchToBody. Set_Vector33 (all,al2, al3, a21, a22,a23, a31, a32, a33);
	return LaunchToBody;
}

//弹体系到发射系，航天一院所提供
C_Matrix33 getMatrix_BodyToLaunch(double IN_dPhi, double IN_dPsi, double IN_dGamma)
{
	C_Matrix33 LaunchToBody;
	C_Matrix33 BodyToLaunch;
	double cPhi = cos(IN_dPhi);
	double sPhi = sin(IN_dPhi);
	double cPsi = cos(IN_dPsi);
	double sPsi = sin(IN_dPsi);
	double cGamma = cos(IN_dGamma);
	double sGamma = sin(IN_dGamma);
	LaunchToBody.a11 = cPhi * cPsi;
	LaunchToBody.a12 = sPhi * cPsi;
	LaunchToBody.a13 = -sPsi;
	LaunchToBody.a21 = cPhi * sPsi * sGamma - sPhi * cGamma;
	LaunchToBody.a22 = sPhi * sPsi * sGamma + cPhi * cGamma;
	LaunchToBody.a23 = cPsi * sGamma;
	LaunchToBody.a31 = cPhi * sPsi * cGamma + sPhi * sGamma;
	LaunchToBody.a32 = sPhi * sPsi * cGamma - cPhi * sGamma;
	LaunchToBody.a33 = cPsi * cGamma;

	BodyToLaunch = BodyToLaunch.fn_CalTran(LaunchToBody);

	return BodyToLaunch;
}

//发射系到速度系，航天一院所提供
C_Matrix33 getMatrix_LaunchToVelocity(double IN_dTheta, double IN_dSegma, double IN_dNu)
{
	C_Matrix33 LaunchToVelocity;
	double cTheta = cos(IN_dTheta);
	double sTheta = sin(IN_dTheta);
	double cSegma = cos(IN_dSegma);
	double sSegma = sin(IN_dSegma);
	double cNu = cos(IN_dNu);
	double sNu = sin(IN_dNu);
	LaunchToVelocity.a11 = cTheta * cSegma;
	LaunchToVelocity.a12 = sTheta * cSegma;
	LaunchToVelocity.a13 = -sSegma;
	LaunchToVelocity.a21 = cTheta * sSegma * sNu - sTheta * cNu;
	LaunchToVelocity.a22 = sTheta * sSegma * sNu + cTheta * cNu;
	LaunchToVelocity.a23 = cSegma * sNu;
	LaunchToVelocity.a31 = cTheta * sSegma * cNu + sTheta * sNu;
	LaunchToVelocity.a32 = sTheta * sSegma * cNu - cTheta * sNu;
	LaunchToVelocity.a33 = cSegma * cNu;

	return LaunchToVelocity;
}

//速度系到发射系，航天一院所提供
C_Matrix33 getMatrix_VelocityToLaunch(double IN_dTheta, double IN_dSegma, double IN_dNu)
{
	C_Matrix33 LaunchToVelocity;
	C_Matrix33 VelocityToLaunch;
	double cTheta = cos(IN_dTheta);
	double sTheta = sin(IN_dTheta);
	double cSegma = cos(IN_dSegma);
	double sSegma = sin(IN_dSegma);
	double cNu = cos(IN_dNu);
	double sNu = sin(IN_dNu);
	LaunchToVelocity.a11 = cTheta * cSegma;
	LaunchToVelocity.a12 = sTheta * cSegma;
	LaunchToVelocity.a13 = -sSegma;
	LaunchToVelocity.a21 = cTheta * sSegma * sNu - sTheta * cNu;
	LaunchToVelocity.a22 = sTheta * sSegma * sNu + cTheta * cNu;
	LaunchToVelocity.a23 = cSegma * sNu;
	LaunchToVelocity.a31 = cTheta * sSegma * cNu + sTheta * sNu;
	LaunchToVelocity.a32 = sTheta * sSegma * cNu - cTheta * sNu;
	LaunchToVelocity.a33 = cSegma * cNu;

	VelocityToLaunch = VelocityToLaunch.fn_CalInv(LaunchToVelocity);
	return VelocityToLaunch;

}

//速度系到弹体系，航天一院所提供
C_Matrix33 getMatrix_VelocityToBody(double IN_dAlpha, double IN_dBeta)
{
	C_Matrix33 VelocityToBody;
	double cAlpha = cos(IN_dAlpha);
	double sAlpha = sin(IN_dAlpha);
	double cBeta = cos(IN_dBeta);
	double sBeta = sin(IN_dBeta);

	VelocityToBody.a11 = cBeta * cAlpha;
	VelocityToBody.a12 = sAlpha;
	VelocityToBody.a13 = -cAlpha * sBeta;
	VelocityToBody.a21 = -cBeta * sAlpha;
	VelocityToBody.a22 = cAlpha;
	VelocityToBody.a23 = sBeta * sAlpha;
	VelocityToBody.a31 = sBeta;
	VelocityToBody.a32 = 0;
	VelocityToBody.a33 = cBeta;

	return VelocityToBody;
}

//弹体系到速度系，航天一院所提供
C_Matrix33 getMatrix_BodyToVelocity(double IN_dAlpha, double IN_dBeta)
{
	C_Matrix33 BodyToVelocity;
	C_Matrix33 VelocityToBody;
	double cAlpha = cos(IN_dAlpha);
	double sAlpha = sin(IN_dAlpha);
	double cBeta = cos(IN_dBeta);
	double sBeta = sin(IN_dBeta);

	VelocityToBody.a11 = cos(IN_dBeta) * cAlpha;
	VelocityToBody.a12 = sAlpha;
	VelocityToBody.a13 = -sin(IN_dBeta) * cAlpha;
	VelocityToBody.a21 = -cos(IN_dBeta) * sAlpha;
	VelocityToBody.a22 = cAlpha;
	VelocityToBody.a23 = sin(IN_dBeta) * sAlpha;
	VelocityToBody.a31 = sin(IN_dBeta);
	VelocityToBody.a32 = 0;
	VelocityToBody.a33 = cos(IN_dBeta);

	BodyToVelocity = BodyToVelocity.fn_CalInv(VelocityToBody);

	return BodyToVelocity;
};

//地心系到发射系，航天一院所提供（未用到）
C_Matrix33 getMatrix_EarthToLaunch(double IN_dLauLon, double IN_dLauLat, double IN_dLauAzi)
{
	C_Matrix33 EarthToLaunch;
	double cLon = cos(IN_dLauLon);
	double sLon = sin(IN_dLauLon);
	double cLat = cos(IN_dLauLat);
	double sLat = sin(IN_dLauLat);
	double cA0 = cos(IN_dLauAzi);
	double sA0 = sin(IN_dLauAzi);
	EarthToLaunch.a11 = -sA0 * sLon - cA0 * sLat * cLon;
	EarthToLaunch.a12 = sA0 * cLon - cA0 * sLat * sLon;
	EarthToLaunch.a13 = cA0 * cLat;
	EarthToLaunch.a21 = cLat * cLon;
	EarthToLaunch.a22 = cLat * sLon;
	EarthToLaunch.a23 = sLat;
	EarthToLaunch.a31 = -cA0 * sLon + sA0 * sLat * cLon;
	EarthToLaunch.a32 = cA0 * cLon + sA0 * sLat * sLon;
	EarthToLaunch.a33 = -sA0 * cLat;

	return EarthToLaunch;
}

//发射系到地心系，航天一院所提供（未用到）
C_Matrix33 getMatrix_LaunchToEarth(double IN_dLauLon, double IN_dLauLat, double IN_dLauAzi)
{
	C_Matrix33 LaunchToEarth;
	C_Matrix33 EarthToLaunch;
	double cLon = cos(IN_dLauLon);
	double sLon = sin(IN_dLauLon);
	double cLat = cos(IN_dLauLat);
	double sLat = sin(IN_dLauLat);
	double cA0 = cos(IN_dLauAzi);
	double sA0 = sin(IN_dLauAzi);
	EarthToLaunch.a11 = -sA0 * sLon - cA0 * sLat * cLon;
	EarthToLaunch.a12 = sA0 * cLon - cA0 * sLat * sLon;
	EarthToLaunch.a13 = cA0 * cLat;
	EarthToLaunch.a21 = cLat * cLon;
	EarthToLaunch.a22 = cLat * sLon;
	EarthToLaunch.a23 = sLat;
	EarthToLaunch.a31 = -cA0 * sLon + sA0 * sLat * cLon;
	EarthToLaunch.a32 = cA0 * cLon + sA0 * sLat * sLon;
	EarthToLaunch.a33 = -sA0 * cLat;

	LaunchToEarth = LaunchToEarth.fn_CalInv(EarthToLaunch);
	return LaunchToEarth;
}

//惯性系到发射系，考虑地球自转
C_Matrix33 getMatrix_InertiaToLaunch(C_Vector3 IN_vLauLLA, double IN_dLauAzi, double IN_dFlyTime)
{
	C_Matrix33 InertiaToLaunch;
	double cLat = cos(IN_vLauLLA.y * DEG2RAD);
	double sLat = sin(IN_vLauLLA.y * DEG2RAD);
	double cA0 = cos(IN_dLauAzi * DEG2RAD);
	double sA0 = sin(IN_dLauAzi * DEG2RAD);
	double cwet = cos(EARTH_WE * IN_dFlyTime);
	double swet = sin(EARTH_WE * IN_dFlyTime);

	InertiaToLaunch.a11 = cA0 * cA0 * cLat * cLat * (1 - cwet) + cwet;
	InertiaToLaunch.a12 = cA0 * sLat * cLat * (1 - cwet) - sA0 * cLat * swet;
	InertiaToLaunch.a13 = -sA0 * cA0 * cLat * cLat * (1 - cwet) - sLat * swet;
	InertiaToLaunch.a21 = cA0 * sLat * cLat * (1 - cwet) + sA0 * cLat * swet;
	InertiaToLaunch.a22 = sLat * sLat * (1 - cwet) + cwet;
	InertiaToLaunch.a23 = -sA0 * sLat * cLat * (1 - cwet) + cA0 * cLat * swet;
	InertiaToLaunch.a31 = -sA0 * cA0 * cLat * cLat * (1 - cwet) + sLat * swet;
	InertiaToLaunch.a32 = -sA0 * sLat * cLat * (1 - cwet) - cA0 * cLat * swet;
	InertiaToLaunch.a33 = sA0 * sA0 * cLat * cLat * (1 - cwet) + cwet;

	return InertiaToLaunch;
}

//发射系到惯性系，考虑地球自转
C_Matrix33 getMatrix_LaunchToInertia(C_Vector3 IN_vLauLLA, double IN_dLauAzi, double IN_dFlyTime)
{
	C_Matrix33 InertiaToLaunch;
	C_Matrix33 LaunchToInertia;
	double cLat = cos(IN_vLauLLA.y * DEG2RAD);
	double sLat = sin(IN_vLauLLA.y * DEG2RAD);
	double cA0 = cos(IN_dLauAzi * DEG2RAD);
	double sA0 = sin(IN_dLauAzi * DEG2RAD);
	double cwet = cos(EARTH_WE * IN_dFlyTime);
	double swet = sin(EARTH_WE * IN_dFlyTime);

	InertiaToLaunch.a11 = cA0 * cA0 * cLat * cLat * (1 - cwet) + cwet;
	InertiaToLaunch.a12 = cA0 * sLat * cLat * (1 - cwet) - sA0 * cLat * swet;
	InertiaToLaunch.a13 = -sA0 * cA0 * cLat * cLat * (1 - cwet) - sLat * swet;
	InertiaToLaunch.a21 = cA0 * sLat * cLat * (1 - cwet) + sA0 * cLat * swet;
	InertiaToLaunch.a22 = sLat * sLat * (1 - cwet) + cwet;
	InertiaToLaunch.a23 = -sA0 * sLat * cLat * (1 - cwet) + cA0 * cLat * swet;
	InertiaToLaunch.a31 = -sA0 * cA0 * cLat * cLat * (1 - cwet) + sLat * swet;
	InertiaToLaunch.a32 = -sA0 * sLat * cLat * (1 - cwet) - cA0 * cLat * swet;
	InertiaToLaunch.a33 = sA0 * sA0 * cLat * cLat * (1 - cwet) + cwet;

	LaunchToInertia = LaunchToInertia.fn_CalInv(InertiaToLaunch);

	return LaunchToInertia;
}

//发射系到弹体系，根据《导弹飞行动力学》所编， 俯仰、偏航、滚转
C_Matrix33 getMatrix_Launch2Body(double phi, double psi, double gamma)
{
	C_Matrix33 Inertia2Body;
	double cPhi = cos(phi);
	double sPhi = sin(phi);
	double cPsi = cos(psi);
	double sPsi = sin(psi);
	double cGamma = cos(gamma);
	double sGamma = sin(gamma);
	Inertia2Body.a11 = cPhi * cPsi;
	Inertia2Body.a12 = sPhi;
	Inertia2Body.a13 = -cPhi * sPsi;
	Inertia2Body.a21 = -sPhi * cPsi * cGamma + sPsi * sGamma;
	Inertia2Body.a22 = cPhi * cGamma;
	Inertia2Body.a23 = sPhi * sPsi * cGamma + cPsi * sGamma;
	Inertia2Body.a31 = sPhi * cPsi * sGamma + sPsi * cGamma;
	Inertia2Body.a32 = -cPhi * sGamma;
	Inertia2Body.a33 = -sPhi * sPsi * sGamma + cPsi * cGamma;

	return Inertia2Body;
}
	
//弹体系到发射系，根据《导弹飞行动力学》所编， 俯仰、偏航、滚转
C_Matrix33 getMatrix_Body2Launch(double phi, double psi, double gamma)
{
	C_Matrix33 Inertia2Body;
	C_Matrix33 Body2Inertia;

	double cPhi = cos(phi);
	double sPhi = sin(phi);
	double cPsi = cos(psi);
	double sPsi = sin(psi);
	double cGamma = cos(gamma);
	double sGamma = sin(gamma);
	Inertia2Body.a11 = cPhi * cPsi;
	Inertia2Body.a12 = sPhi;
	Inertia2Body.a13 = -cPhi * sPsi;
	Inertia2Body.a21 = -sPhi * cPsi * cGamma + sPsi * sGamma;
	Inertia2Body.a22 = cPhi * cGamma;
	Inertia2Body.a23 = sPhi * sPsi * cGamma + cPsi * sGamma;
	Inertia2Body.a31 = sPhi * cPsi * sGamma + sPsi * cGamma;
	Inertia2Body.a32 = -cPhi * sGamma;
	Inertia2Body.a33 = -sPhi * sPsi * sGamma + cPsi * cGamma;

	Body2Inertia = Body2Inertia.fn_CalInv(Inertia2Body);

	return Body2Inertia;
}

//弹道系到速度系，根据《导弹飞行动力学》所编，倾侧角
C_Matrix33 getMatrix_Dandao2Velocity(double nu)
{
	C_Matrix33 dandao2velocity;

	dandao2velocity.a11 = 1;
	dandao2velocity.a12 = 0;
	dandao2velocity.a13 = 0;
	dandao2velocity.a21 = 0;
	dandao2velocity.a22 = cos(nu);
	dandao2velocity.a23 = sin(nu);
	dandao2velocity.a31 = 0;
	dandao2velocity.a32 = -sin(nu);
	dandao2velocity.a33 = cos(nu);

	return dandao2velocity;
}

//速度系到弹道系，根据《导弹飞行动力学》所编，倾侧角
C_Matrix33 getMatrix_Velocity2Dandao(double nu)
{
	C_Matrix33 dandao2velocity;
	C_Matrix33 velocity2dandao;

	dandao2velocity.a11 = 1;
	dandao2velocity.a12 = 0;
	dandao2velocity.a13 = 0;
	dandao2velocity.a21 = 0;
	dandao2velocity.a22 = cos(nu);
	dandao2velocity.a23 = sin(nu);
	dandao2velocity.a31 = 0;
	dandao2velocity.a32 = -sin(nu);
	dandao2velocity.a33 = cos(nu);

	velocity2dandao = velocity2dandao.fn_CalInv(dandao2velocity);

	return velocity2dandao;
}

//速度系到弹体系，根据《导弹飞行动力学》所编，攻角，侧滑角
C_Matrix33 getMatrix_Velocity2Body(double IN_dAlpha, double IN_dBeta)
{
	C_Matrix33 VelocityToBody;
	double cAlpha = cos(IN_dAlpha);
	double sAlpha = sin(IN_dAlpha);
	double cBeta = cos(IN_dBeta);
	double sBeta = sin(IN_dBeta);

	VelocityToBody.a11 = cBeta * cAlpha;
	VelocityToBody.a12 = sAlpha;
	VelocityToBody.a13 = -cAlpha * sBeta;
	VelocityToBody.a21 = -cBeta * sAlpha;
	VelocityToBody.a22 = cAlpha;
	VelocityToBody.a23 = sBeta * sAlpha;
	VelocityToBody.a31 = sBeta;
	VelocityToBody.a32 = 0;
	VelocityToBody.a33 = cBeta;

	return VelocityToBody;
}

//弹体系到速度系，根据《导弹飞行动力学》所编，攻角，侧滑角
C_Matrix33 getMatrix_Body2Velocity(double IN_dAlpha, double IN_dBeta)
{
	C_Matrix33 BodyToVelocity;
	C_Matrix33 VelocityToBody;
	double cAlpha = cos(IN_dAlpha);
	double sAlpha = sin(IN_dAlpha);
	double cBeta = cos(IN_dBeta);
	double sBeta = sin(IN_dBeta);

	VelocityToBody.a11 = cos(IN_dBeta) * cAlpha;
	VelocityToBody.a12 = sAlpha;
	VelocityToBody.a13 = -sin(IN_dBeta) * cAlpha;
	VelocityToBody.a21 = -cos(IN_dBeta) * sAlpha;
	VelocityToBody.a22 = cAlpha;
	VelocityToBody.a23 = sin(IN_dBeta) * sAlpha;
	VelocityToBody.a31 = sin(IN_dBeta);
	VelocityToBody.a32 = 0;
	VelocityToBody.a33 = cos(IN_dBeta);

	BodyToVelocity = BodyToVelocity.fn_CalInv(VelocityToBody);

	return BodyToVelocity;
}

//ECEF坐标系转东北天，输入的LLA单位是度
C_Vector3 ECEF2ENU(C_Vector3 in_Launchpoint_ECEF, C_Vector3 in_LLA, C_Vector3 in_ECEF)
{
	C_Vector3 dEcef = in_ECEF - in_Launchpoint_ECEF; //平移
	double slon = sin(in_LLA.x * DEG2RAD);
	double clon = cos(in_LLA.x * DEG2RAD);
	double slat = sin(in_LLA.y * DEG2RAD);
	double clat = cos(in_LLA.y * DEG2RAD);
	C_Matrix33 ECEFtoENU;
	ECEFtoENU.a11 = -slon;
	ECEFtoENU.a12 = clon;
	ECEFtoENU.a13 = 0;
	ECEFtoENU.a21 = -slat * clon;
	ECEFtoENU.a22 = -slat * slon;
	ECEFtoENU.a23 = clat;
	ECEFtoENU.a31 = clat * clon;
	ECEFtoENU.a32 = clat * slon;
	ECEFtoENU.a33 = slat;

	return ECEFtoENU * dEcef; //旋转

}

//东北天到发射系（苏联坐标系）
C_Vector3 ENU2Launch(double in_Launchalpha, C_Vector3 in_ENU)
{
	C_Vector3 launch;
	C_Matrix33 ENUto2launch;
	ENUto2launch.a11 = -sin(in_Launchalpha);
	ENUto2launch.a12 = cos(in_Launchalpha);
	ENUto2launch.a13 = 0;
	ENUto2launch.a21 = 0;
	ENUto2launch.a22 = 0;
	ENUto2launch.a23 = 1;
	ENUto2launch.a31 = cos(in_Launchalpha);
	ENUto2launch.a32 = sin(in_Launchalpha);
	ENUto2launch.a33 = 0;

	launch = ENUto2launch * in_ENU;
	return launch;

}

//发射系转东北天，输入的LLA单位是度
C_Vector3 Launch2ENU(double in_Launchalpha, C_Vector3 in_launch)
{
	//东北天转为北天东
	C_Vector3 enu;
	C_Matrix33 ENUto2launch;
	C_Matrix33 launchto2ENU;
	ENUto2launch.a11 = -sin(in_Launchalpha);
	ENUto2launch.a12 = cos(in_Launchalpha);
	ENUto2launch.a13 = 0;
	ENUto2launch.a21 = 0;
	ENUto2launch.a22 = 0;
	ENUto2launch.a23 = 1;
	ENUto2launch.a31 = cos(in_Launchalpha);
	ENUto2launch.a32 = sin(in_Launchalpha);
	ENUto2launch.a33 = 0;

	launchto2ENU = launchto2ENU.fn_CalTran(ENUto2launch);

	enu = launchto2ENU * in_launch;
	return enu;

}

//ECEF转发射系，输入的LLA单位是度
C_Vector3 ECEF2Launch(C_Vector3 in_Launchpoint_ECEF, C_Vector3 in_LLA, C_Vector3 in_ECEF, double in_Launchalpha)
{
	C_Vector3 enu = ECEF2ENU(in_Launchpoint_ECEF, in_LLA, in_ECEF);
	C_Vector3 launch = ENU2Launch(in_Launchalpha, enu);
	return launch;

}

//发射系转ECEF，输入的LLA单位是度
C_Vector3 Launch2ECEF(C_Vector3 in_Launchpoint_ECEF, C_Vector3 in_LLA, C_Vector3 in_launch, double in_Launchalpha)
{
	C_Vector3 enu = Launch2ENU(in_Launchalpha, in_launch);
	C_Vector3 ecef = ENU2ECEF(in_Launchpoint_ECEF, in_LLA, enu);

	return ecef;
}

//东北天转ECEF，输入的LLA单位是度
C_Vector3 ENU2ECEF(C_Vector3 in_Launchpoint_ECEF, C_Vector3 in_LLA, C_Vector3 in_ENU)
{
	C_Vector3 ECEF;
	double slon = sin(in_LLA.x * DEG2RAD);
	double clon = cos(in_LLA.x * DEG2RAD);
	double slat = sin(in_LLA.y * DEG2RAD);
	double clat = cos(in_LLA.y * DEG2RAD);
	C_Matrix33 ECEFtoENU;
	C_Matrix33 ENUtoECEF;
	ECEFtoENU.a11 = -slon;
	ECEFtoENU.a12 = clon;
	ECEFtoENU.a13 = 0;
	ECEFtoENU.a21 = -slat * clon;
	ECEFtoENU.a22 = -slat * slon;
	ECEFtoENU.a23 = clat;
	ECEFtoENU.a31 = clat * clon;
	ECEFtoENU.a32 = clat * slon;
	ECEFtoENU.a33 = slat;

	ENUtoECEF = ENUtoECEF.fn_CalTran(ECEFtoENU);
	ECEF = ENUtoECEF * in_ENU;  //旋转

	return ECEF + in_Launchpoint_ECEF; //平移

};

//绕x轴旋转矩阵
C_Matrix33 getMatrix_RotateX(double IN_dAngle)
{
	C_Matrix33 mRotateX;
	mRotateX.a11 = 1;
	mRotateX.a12 = 0;
	mRotateX.a13 = 0;
	mRotateX.a21 = 0;
	mRotateX.a22 = cos(IN_dAngle);
	mRotateX.a23 = sin(IN_dAngle);
	mRotateX.a31 = 0;
	mRotateX.a32 = -sin(IN_dAngle);
	mRotateX.a33 = cos(IN_dAngle);

	return mRotateX;
}

//绕y轴旋转矩阵
C_Matrix33 getMatrix_RotateY(double IN_dAngle)
{
	C_Matrix33 mRotateY;
	mRotateY.a11 = cos(IN_dAngle);
	mRotateY.a12 = 0;
	mRotateY.a13 = -sin(IN_dAngle);
	mRotateY.a21 = 0;
	mRotateY.a22 = 1;
	mRotateY.a23 = 0;
	mRotateY.a31 = sin(IN_dAngle);
	mRotateY.a32 = 0;
	mRotateY.a33 = cos(IN_dAngle);

	return mRotateY;
}

//绕z轴旋转矩阵
C_Matrix33 getMatrix_RotateZ(double IN_dAngle)
{
	C_Matrix33 mRotateZ;
	mRotateZ.a11 = cos(IN_dAngle);
	mRotateZ.a12 = sin(IN_dAngle);
	mRotateZ.a13 = 0;
	mRotateZ.a21 = -sin(IN_dAngle);
	mRotateZ.a22 = cos(IN_dAngle);
	mRotateZ.a23 = 0;
	mRotateZ.a31 = 0;
	mRotateZ.a32 = 0;
	mRotateZ.a33 = 1;

	return mRotateZ;
}

//ECI到ECEF坐标系
C_Matrix33 getMatrix_ECItoECEF(double IN_dFlyTime)
{
	C_Matrix33 InertiaToEarth;
	double wet = EARTH_WE * IN_dFlyTime;
	InertiaToEarth = getMatrix_RotateZ(wet);
	return InertiaToEarth;
}

//ECI到ECEF坐标系
C_Matrix33 getMatrix_ECEFtoECI(double IN_dFlyTime)
{
	C_Matrix33 EarthToInertia;
	EarthToInertia = EarthToInertia.fn_CalInv(getMatrix_ECItoECEF(IN_dFlyTime));
	return EarthToInertia;
}

//ECI到发射系
C_Vector3 Transform_ECItoLaunch(C_Vector3 IN_vLauLLA, double IN_dLauAzi, double IN_dFlyTime, C_Vector3 IN_vPosECI)
{
	C_Vector3 vPosECEF = getMatrix_ECItoECEF(IN_dFlyTime) * IN_vPosECI;
	C_Matrix33 mECEFtoLau = getMatrix_EarthToLaunch(IN_vLauLLA.x * DEG2RAD, IN_vLauLLA.y * DEG2RAD, IN_dLauAzi * DEG2RAD);
	C_Vector3 vPosLauPoint((EARTH_RAVG + IN_vLauLLA.z) * cos(IN_vLauLLA.x * DEG2RAD) * cos(IN_vLauLLA.y * DEG2RAD),
		(EARTH_RAVG + IN_vLauLLA.z) * sin(IN_vLauLLA.x * DEG2RAD) * cos(IN_vLauLLA.y * DEG2RAD),
		(EARTH_RAVG + IN_vLauLLA.z) * sin(IN_vLauLLA.y * DEG2RAD));
	C_Vector3 vPosLau = mECEFtoLau * (vPosECEF - vPosLauPoint);
	return vPosLau;
}

//84LLA经纬高到ECEF坐标系,输入是度
C_Vector3 Transform_LLA2ECEF(C_Vector3 vLLA)
{
	C_Vector3 vPosECEF;
	double lon, lat, h;
	lon = vLLA.x * DEG2RAD;
	lat = vLLA.y * DEG2RAD;
	h = vLLA.z;
	vPosECEF.x = (EARTH_RAVG + h) * cos(lon) * cos(lat);
	vPosECEF.y = (EARTH_RAVG + h) * sin(lon) * cos(lat);
	vPosECEF.z = (EARTH_RAVG + h) * sin(lat);
	return vPosECEF;
}

//ECEF坐标系到84LLA经纬高，输出是弧度
C_Vector3 Transform_ECEF2LLA(C_Vector3 IN_vEarthPose)
{
	C_Vector3 vEarthLLA;
	//计算地心系下xY方向位置平方和
	double dEarthXYabs = sqrt((IN_vEarthPose.x * IN_vEarthPose.x) + (IN_vEarthPose.y * IN_vEarthPose.y));
	//计算当前地心纬度
	vEarthLLA.y = atan2(IN_vEarthPose.z, dEarthXYabs);
	//分情况计算经度
	if (IN_vEarthPose.x > 0 && IN_vEarthPose.y >= 0)
	{
		vEarthLLA.x = acos(IN_vEarthPose.x / dEarthXYabs);
	}
	if (IN_vEarthPose.x <= 0 && IN_vEarthPose.y > 0)
	{
		vEarthLLA.x = PI - acos(fabs(IN_vEarthPose.x) / dEarthXYabs);
	}
	if (IN_vEarthPose.x < 0 && IN_vEarthPose.y <= 0)
		vEarthLLA.x = -PI + acos(fabs(IN_vEarthPose.x) / dEarthXYabs);
	if (IN_vEarthPose.x >= 0 && IN_vEarthPose.y < 0)
		vEarthLLA.x = -acos(fabs(IN_vEarthPose.x) / dEarthXYabs);
	if (IN_vEarthPose.x == 0 && IN_vEarthPose.y == 0)
		vEarthLLA.x = 0;

	//地心系下当前弹下点地球半径
	double dEarthR = EARTH_RAVG;
	//计算高度
	double dEarthHeight = IN_vEarthPose.Mod() - dEarthR;
	//高度输出
	vEarthLLA.z = dEarthHeight;
	return vEarthLLA;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//三维数组求差
void minusVec2(double Vec2Input1[2], double Vec2Input2[2], double Vec2Output[2])
{
	for (int i = 0; i < 2; i++)
	{
		Vec2Output[i] = Vec2Input1[i] - Vec2Input2[i];
	}
}

//计算地球上两点的球面距离
double g_fnCalSphericalDist(double d_lon1, double d_lat1, double d_lon2, double d_lat2)
{
	double d_sphericaldist;
	d_sphericaldist = EARTH_RAVG * acos(cos(d_lat1 * DEG2RAD) * cos(d_lat2 * DEG2RAD) * cos(d_lon1 * DEG2RAD - d_lon2 * DEG2RAD) +
		sin(d_lat1 * DEG2RAD) * sin(d_lat2 * DEG2RAD));
	return d_sphericaldist;
}

//验证double类型的两个数是否相等
bool g_fnCheckDoubleEqual(double xl, double x2)
{
	if (fabs(xl - x2) < SMALLNUM)
	{
		return true;
	}
	else
	{
		return false;
	}
}

double g_fnCalSphericalAzimuth(double d_laulon, double d_laulat, double d_tarlon, double d_tarlat)
{
	double d_laulon_rad = d_laulon * DEG2RAD;
	double d_laulat_rad = d_laulat * DEG2RAD;
	double d_tarlon_rad = d_tarlon * DEG2RAD;
	double d_tarlat_rad = d_tarlat * DEG2RAD;

	double m_d_A0 = sin(d_laulat_rad) * sin(d_tarlat_rad) + cos(d_laulat_rad) * cos(d_tarlat_rad) * cos(d_tarlon_rad - d_laulon_rad);
	m_d_A0 = sqrt(1 - m_d_A0 * m_d_A0);
	m_d_A0 = cos(d_tarlat_rad) * sin(d_tarlon_rad - d_laulon_rad) / m_d_A0;
	m_d_A0 = asin(m_d_A0) * RAD2DEG;
	if (_finite(m_d_A0) != 1)
	{
		if (d_laulon_rad < d_tarlon_rad)
		{
			m_d_A0 = 90.0;
		}
		else
		{
			m_d_A0 = 270.0;
		}
	}
	double d_dLon_rad = d_tarlon_rad - d_laulon_rad;
	double d_dLat_rad = d_tarlat_rad - d_laulat_rad;

	if (d_dLon_rad > 0 && d_dLat_rad > 0)
	{

	}
	else if (d_dLon_rad < 0 && d_dLat_rad > 0)
	{
		m_d_A0 += 360;
	}

	else if (d_dLon_rad < 0 && d_dLat_rad < 0)
	{
		m_d_A0 = 180 - m_d_A0;
	}
	else if (d_dLon_rad > 0 && d_dLat_rad < 0)
	{
		m_d_A0 = 180 - m_d_A0;
	}

	m_d_A0 *= DEG2RAD;
	return m_d_A0;

}

//根据高度和温度计算大气压强
double g_fnCalAtomspherePressure(double IN_dH, double IN_dT)
{
	const double ATMO_P0 = 101325.0;                                 // 标准大气压(pa)
	const double ATMO_K = 273.15;                                   // 摄氏开氏温度互换
	return  ATMO_P0 / exp(IN_dH / (18400 * (1 + IN_dT / ATMO_K)));
}

//计算一个数组前n个元素的和
double g_fnCalSum(const double* IN_p_arr_dx, int IN_iDatalength)
{
	double d_sum = 0;
	for (int i = 0; i < IN_iDatalength; i++)
	{
		d_sum += IN_p_arr_dx[i];
	}
	return d_sum;
}

//计算一个数组前n个元素的平均值
double g_fnCalAvg(const double* IN_p_arr_dx, int IN_iDatalength)
{
	double d_avg = 0;
	for (int i = 0; i < IN_iDatalength; i++)
	{
		d_avg += IN_p_arr_dx[i];
	}
	return d_avg / IN_iDatalength;
}

//计算一个数组前n个元素的最大值
double g_fnCalMax(const double* IN_p_arr_dArray, int IN_iDatalength)
{
	double dMax = IN_p_arr_dArray[0];
	for (int i = 0; i < IN_iDatalength; i++)
	{
		if (IN_p_arr_dArray[i] > dMax)
		{
			dMax = IN_p_arr_dArray[i];
		}
	}
	return dMax;
}

//计算两个三维向量的夹角
double g_fnCalVectorAngle(C_Vector3 IN_vVec1, C_Vector3 IN_vVec2)
{
	double d_mod1 = IN_vVec1.Mod();
	double d_mod2 = IN_vVec2.Mod();

	double d_cosAgl = (IN_vVec1.x * IN_vVec2.x + IN_vVec1.y * IN_vVec2.y + IN_vVec1.z * IN_vVec2.z) / (d_mod1 * d_mod2);

	double d_Agl = acos(d_cosAgl);

	return d_Agl;

}

//double类型的数向上取整
int g_fnRoundDouble(double number)
{
	return (int)((number > 0.0) ? (number + 0.5) : (number - 0.5));
}

//计算LLA坐标系下的两点距离
double g_fnCalLLADist(C_Vector3 IN_vLLA1, C_Vector3 IN_vLLA2)
{
	C_Vector3 v_earth1;//点1地心系位置
	double d_centredist1 = IN_vLLA1.z + EARTH_RAVG;//点1地心距
	v_earth1.z = d_centredist1 * sin(IN_vLLA1.y * DEG2RAD);
	v_earth1.x = d_centredist1 * cos(IN_vLLA1.y * DEG2RAD) * cos(IN_vLLA1.x * DEG2RAD);
	v_earth1.y = d_centredist1 * cos(IN_vLLA1.y * DEG2RAD) * sin(IN_vLLA1.x * DEG2RAD);

	C_Vector3 v_earth2;//点1地心系位置
	double d_centredist2 = IN_vLLA2.z + EARTH_RAVG;//点1地心距
	v_earth2.z = d_centredist2 * sin(IN_vLLA2.y * DEG2RAD);
	v_earth2.x = d_centredist2 * cos(IN_vLLA2.y * DEG2RAD) * cos(IN_vLLA2.x * DEG2RAD);
	v_earth2.y = d_centredist2 * cos(IN_vLLA2.y * DEG2RAD) * sin(IN_vLLA2.x * DEG2RAD);

	return (v_earth1 - v_earth2).Mod();

}

//简单的一维直线插值，计算dx（ddelta）的y值
double g_fnOneDimInterp(double IN_dLBoundary, double IN_dRBoundary, double IN_dStep, double IN_dDelta)
{
	double dK = (IN_dRBoundary - IN_dLBoundary) / (IN_dStep + EPS);
	double dOut = dK * IN_dDelta + IN_dLBoundary;

	return dOut;
}

void g_fnThreeDimInterp_Parabolic(const double IN_p_dTXYZ0[4], const double IN_p_dTXYZ1[4], const double IN_p_dTXYZ2[4], double IN_dInterpVal, double OUT_p_dTXYZ[4], double OUT_p_dTDXDYDZ[4])
{
	double dA = 0;
	double dB = 0;
	double dC = 0;

	OUT_p_dTXYZ[0] = IN_dInterpVal;
	OUT_p_dTDXDYDZ[0] = IN_dInterpVal;

	for (int i = 1; i < 4; i++)
	{
		dA = (IN_p_dTXYZ2[i] -
			IN_p_dTXYZ2[0] * (IN_p_dTXYZ0[i] - IN_p_dTXYZ1[i]) / (IN_p_dTXYZ0[0] - IN_p_dTXYZ1[0]) -
			IN_p_dTXYZ0[i] +
			IN_p_dTXYZ0[0] * (IN_p_dTXYZ0[i] - IN_p_dTXYZ1[i]) / (IN_p_dTXYZ0[0] - IN_p_dTXYZ1[0])) /
			(IN_p_dTXYZ2[0] * IN_p_dTXYZ2[0] -
				IN_p_dTXYZ2[0] * (IN_p_dTXYZ0[0] + IN_p_dTXYZ1[0]) -
				IN_p_dTXYZ0[0] * IN_p_dTXYZ0[0] +
				IN_p_dTXYZ0[0] * (IN_p_dTXYZ0[0] + IN_p_dTXYZ1[0]));
		dB = (IN_p_dTXYZ0[i] - IN_p_dTXYZ1[i]) / (IN_p_dTXYZ0[0] - IN_p_dTXYZ1[0]) - (IN_p_dTXYZ0[0] + IN_p_dTXYZ1[0]) * dA;
		dC = IN_p_dTXYZ0[i] - dA * IN_p_dTXYZ0[0] * IN_p_dTXYZ0[0] - dB * IN_p_dTXYZ0[0];

		OUT_p_dTXYZ[i] = dA * IN_dInterpVal * IN_dInterpVal + dB * IN_dInterpVal + dC;
		OUT_p_dTDXDYDZ[i] = 2 * dA * IN_dInterpVal + dB;
	}
	return;

}

void g_fnOrdinaryLeastSquares_Polynomial(int IN_iInputLength, const vector <double>& IN_p_arr_dX, const vector<double >& IN_p_arr_dY, int IN_iOrder,int IN_iOutputLength, const vector <double>& IN_p_arr_dX_esti, vector <double>& OUT_p_arr_dY_esti)
{
		vector <vector <double>> mat_A;
		vector <vector <double>> mat_AT;//A的转置

		//求解矩阵A
		//矩阵的行，行数为输入项数量
		for (int i = 0; i < IN_iInputLength; i++)
		{
			//矩阵每一行的各个列
			vector <double> singleline;
			for (int j = 0; j <= IN_iOrder; j++)
			{
				double dA_ij = 1;
				for (int k = 0; k < j; k++)
				{
					dA_ij *= IN_p_arr_dX[i];
				}
				singleline.push_back(dA_ij);
			}

			mat_A.push_back(singleline);
			singleline.clear();
		}

		//求解矩阵A的转置
		//转置矩阵的行，行数为A的列数
		for (int i = 0; i <= IN_iOrder; i++)
		{
			vector <double> singleline;
			for (int j = 0; j < IN_iInputLength; j++)
			{
				double dAT_ij = mat_A[j][i];
				singleline.push_back(dAT_ij);
			}

			mat_AT.push_back(singleline);
			singleline.clear();
		}

		//求解AT*A
		vector <vector <double>> mat_ATA;
		g_fnCalMatrixMulti(mat_AT, mat_A, mat_ATA);

		//求解AT*b
		vector <vector <double>> mat_b;
		for (int i = 0; i < IN_iInputLength; i++)
		{
			vector <double> singleline;
			singleline.push_back(IN_p_arr_dY[i]);
			mat_b.push_back(singleline);
			singleline.clear();
		}

		vector <vector <double>> mat_ATb;
		g_fnCalMatrixMulti(mat_AT, mat_b, mat_ATb);

		//求解Inv(AT*A)
		vector <vector <double>> mat_InvATA;
		g_fnCalMatrixInv(mat_ATA, mat_InvATA);

		//求解多项式系数矩阵Iv(AT*A)*b
		vector <vector <double>> mat_InvATA_b;
		g_fnCalMatrixMulti(mat_InvATA, mat_ATb, mat_InvATA_b);

		//估计
		int iRow_InvATA_b = 0;
		int iCol_InvATA_b = 0;
		bool bIfError = false;

		g_fnGetMatrixRowCol(mat_InvATA_b, &iRow_InvATA_b, &iCol_InvATA_b, &bIfError);
		if (mat_InvATA_b.size() != IN_iOrder + 1 || (bIfError || iCol_InvATA_b != 1))
		{
			cout << "ERROR: ESTI COEF NUM NOT EQUAL TO ORDER" << endl;
			mat_A.clear();
			mat_AT.clear();
			mat_ATA.clear();
			mat_ATb.clear();
			mat_b.clear();
			mat_InvATA.clear();
			mat_InvATA_b.clear();
			return;
		}
		OUT_p_arr_dY_esti.clear();
		for (int i = 0; i < IN_iOutputLength; i++)
		{
			double dResult_i = 0;
			for (unsigned int j = 0; j < mat_InvATA_b.size(); j++)
			{
				dResult_i += mat_InvATA_b[j][0] * g_fnIntPow(IN_p_arr_dX_esti[i], j);
			}
			OUT_p_arr_dY_esti.push_back(dResult_i);
		}

		//RELEASE MEMORY
		mat_A.clear();
		mat_AT.clear();
		mat_ATA.clear();
		mat_ATb.clear();
		mat_b.clear();
		mat_InvATA.clear();
		mat_InvATA_b.clear();

		return;
}

/*求 n 阶矩阵的行列式函数*/
double g_fnCalMatrixDet(const vector <vector <double>>& IN_p_vec_da)
{
	vector <vector <double>> a;
	for (unsigned int iA = 0; iA < IN_p_vec_da.size(); iA++)
	{
		vector <double> singleline;
		for (unsigned int j = 0; j < IN_p_vec_da[iA].size(); j++)
		{
			singleline.push_back(IN_p_vec_da[iA][j]);
		}
		a.push_back(singleline);
		singleline.clear();
	}
	//检验行列数是否相等
	bool bIfRowColEqual = true;
	for (unsigned int iA = 0; iA < a.size(); iA++)
	{
		if (a.size() != a[iA].size())
		{
			bIfRowColEqual = false;
			break;
		}
	}
	if (!bIfRowColEqual)
	{
		cout << "ERROR: MATRIX LINE WIDTH NOT EQUAL" << endl;
		return 0.0;
	}

	int n = a.size();
	int i, j, k, m, r, p = 0;
	double t, max, t1;
	for (i = 0; i < n; i++)
	{
		for (int iI = 0; iI < n; iI++)
		{
			a[i].push_back(0.0);
		}
		a[i][n + i] = 1; /*构造增广矩阵[A I]*/
	}
	for (k = 0; k < n - 1; k++)/*列主元高斯消去法，将 A 阵化为上三角矩阵*/
	{
		max = a[k][k];
		r = k;/*选主元*/
		for (m = k; m < n; m++)/*找主元及其所在的行*/
		{
			if (fabs(a[m][k]) > max)
			{
				max = a[m][k];
				r = m;
			}
			if (r != k)/*计算行或列交换的次数*/
			{
				p = p + 1;
			}
			for (j = 0; j < 2 * n; j++)/*第 k 行与第 r 行进行交换*/
			{
				t = a[k][j];
				a[k][j] = a[r][j];
				a[r][j] = t;
			}
			for (i = k + 1; i < n; i++)/*进行消元过程*/
			{
				t = a[i][k] / a[k][k];
				for (j = k; j < 2 * n; j++)
				{
					a[i][j] = a[i][j] - t * a[k][j];
				}

			}

		}

	}
	t1 = 1;/*求(-1)~p 的值*/
	for (i = 0; i < p; i++)
	{
		t1 = t1 * (-1);
	}
	t = 1.0;
	for (i = 0; i < n; i++)/*根据性质求解行列式 det(A)*/
	{
		t = t * a[i][i];
	}

	//printf("det(A)=%f\n,t1*t);/*输出 A 阵的行列式*
	a.clear();
	return (t1 * t);
}

//计算每一行每一列的每个元素所对应的余子式，组成A*
void g_fnCalMatrixAStar(const vector <vector <double>>& IN_p_vec_da, vector<vector <double>>& ans)
{
	vector <vector <double>> arcs;
	for (unsigned int iA = 0; iA < IN_p_vec_da.size(); iA++)
	{
		vector <double> singleline;
		for (unsigned int j = 0; j < IN_p_vec_da[iA].size(); j++)
		{
			singleline.push_back(IN_p_vec_da[iA][j]);
		}
		arcs.push_back(singleline);
		singleline.clear();
	}
	//检验行列数是否相等
	bool bIfRowColEqual = true;
	for (unsigned int iA = 0; iA < arcs.size(); iA++)
	{
		if (arcs.size() != arcs[iA].size())
		{
			bIfRowColEqual = false;
			break;
		}
	}
	if (!bIfRowColEqual)
	{
		cout << "ERROR: MATRIX LINE WIDTH NOT EQUAL" << endl;
		return;
	}

	int n = arcs.size();

	ans.clear();
	for (int iA = 0; iA < n; iA++)
	{
		vector <double> singleline;
		for (int jA = 0; jA < n; jA++)
		{
			singleline.push_back(0.0);
		}
		ans.push_back(singleline);
		singleline.clear();
	}

	if (n == 1)
	{
		ans[0][0] = 1.0;
		return;
	}

	int i, j, k, t;
	//double temp[N] [N] :

	vector <vector <double>> temp;
	for (unsigned int iA = 0; iA < arcs.size(); iA++)
	{
		vector <double> singleline;
		for (unsigned int jA = 0; jA < arcs.size(); jA++)
		{
			singleline.push_back(0.0);
		}
		temp.push_back(singleline);
		singleline.clear();
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < n - 1; k++)
			{
				for (t = 0; t < n - 1; t++)
				{
					temp[k][t] = arcs[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
				}
			}
			vector <vector <double>> tempAij;
			for (int iA = 0; iA < n - 1; iA++)
			{
				vector <double> singleline;
				for (int jA = 0; jA < n - 1; jA++)
				{
					singleline.push_back(temp[iA][jA]);
				}
				tempAij.push_back(singleline);
				singleline.clear();
			}
			ans[j][i] = g_fnCalMatrixDet(tempAij); //此处顺便进行了转置
			tempAij.clear();
			if ((i + j) % 2 == 1)
			{
				ans[j][i] = -ans[j][i];
			}
		}
	}
	arcs.clear();
	temp.clear();
	return;
}

//得到给定矩阵src的逆矩阵保存到des中
bool g_fnCalMatrixInv(const vector <vector <double>>& IN_p_vec_da, vector<vector <double >>& des)
{
	vector <vector <double>> src;
	for (unsigned int iA = 0; iA < IN_p_vec_da.size(); iA++)
	{
		vector <double> singleline;
		for (unsigned int j = 0; j < IN_p_vec_da[iA].size(); j++)
		{
			singleline.push_back(IN_p_vec_da[iA][j]);
		}
		src.push_back(singleline);
		singleline.clear();
	}
	//检验行列数是否相等
	bool bIfRowColEqual = true;
	for (unsigned int iA = 0; iA < src.size(); iA++)
	{
		if (src.size() != src[iA].size())
		{
			bIfRowColEqual = false;
			break;
		}
	}

	if (!bIfRowColEqual)
	{
		cout << "ERROR: MATRIX LINE WIDTH NOT EQUAL" << endl;
		src.clear();
		return false;
	}

	int n = src.size();
	double flag = g_fnCalMatrixDet(src);
	vector <vector <double>> t;
	for (unsigned int iA = 0; iA < IN_p_vec_da.size(); iA++)
	{
		vector <double> singleline;
		for (int j = 0; j < IN_p_vec_da.size(); j++)
		{
			singleline.push_back(0.0);
		}
		t.push_back(singleline);
		singleline.clear();
	}

	if (flag == 0)
	{
		cout << "原矩阵行列式为0，无法求逆。请重新运行" << endl;
		src.clear();
		t.clear();
		return false;//如果算出矩阵的行列式为0，则不往下进行
	}

	else
	{
		g_fnCalMatrixAStar(src, t);

		des.clear();
		for (int iA = 0; iA < n; iA++)
		{
			vector <double> singleline;
			for (int jA = 0; jA < n; jA++)
			{
				singleline.push_back(0.0);
			}

			des.push_back(singleline);
			singleline.clear();
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				des[i][j] = t[i][j] / flag;
			}

		}
	}
	src.clear();
	t.clear();
	return true;

}

//求矩阵的行列数
void g_fnGetMatrixRowCol(const vector <vector <double>>& IN_p_vec_dMatrix, int* OUT_p_iRow, int* OUT_p_iCol, bool* OUT_p_bIfError)
{
	*OUT_p_iRow = 0;
	*OUT_p_iCol = 0;
	*OUT_p_bIfError = false;

	unsigned int uiRow = IN_p_vec_dMatrix.size();
	*OUT_p_iRow = uiRow;
	if (uiRow == 0)
	{
		cout << "MATRIX ROW NUM ERROR" << endl;
		*OUT_p_iRow = 0;
		*OUT_p_iCol = 0;
		*OUT_p_bIfError = true;
		return;
	}
	for (unsigned int i = 0; i < uiRow; i++)
	{
		if (IN_p_vec_dMatrix[i].size() == 0)
		{
			*OUT_p_bIfError = true;
			break;
		}

		if (i == 0)
		{
			*OUT_p_iCol = IN_p_vec_dMatrix[i].size();
		}
		else
		{
			if (*OUT_p_iCol != IN_p_vec_dMatrix[i].size())
			{
				cout << "MATRIX COL NUM ERROR" << endl;
				*OUT_p_bIfError = true;
				break;
			}
		}
	}
	return;
}

//矩阵相乘
void g_fnCalMatrixMulti(const vector<vector <double>>& IN_p_vec_dMat1, const vector<vector<double>>& IN_p_vec_dMat2, vector <vector <double>>& OUT_p_vec_dResult)
{
	int iRow1 = 0;
	int iCol1 = 0;
	bool bIfError1 = false;
	int iRow2 = 0;
	int iCol2 = 0;
	bool bIfError2 = false;

	g_fnGetMatrixRowCol(IN_p_vec_dMat1, &iRow1, &iCol1, &bIfError1);
	g_fnGetMatrixRowCol(IN_p_vec_dMat2, &iRow2, &iCol2, &bIfError2);

	if (bIfError1 || bIfError2)
	{
		OUT_p_vec_dResult.clear();
		return;
	}
	if (iCol1 != iRow2)
	{
		cout << "MATRRIX1 COL NOT EQUAL TO MATRIX2 ROW" << endl;
		OUT_p_vec_dResult.clear();
		return;
	}

	OUT_p_vec_dResult.clear();
	for (int i = 0; i < iRow1; i++)
	{
		vector <double> singleline;
		for (int j = 0; j < iCol2; j++)
		{
			double dResult_ij = 0;
			for (int k = 0; k < iCol1; k++)
			{
				dResult_ij += IN_p_vec_dMat1[i][k] * IN_p_vec_dMat2[k][j];
			}
			singleline.push_back(dResult_ij);
		}
		OUT_p_vec_dResult.push_back(singleline);
		singleline.clear();
	}
	return;
};

//检查输入值是否在界限内
bool g_fnCheckInRange(double IN_dLeftBoundary, double IN_dVar, double IN_dRightBoundary)
{
	if (IN_dVar > IN_dLeftBoundary && IN_dVar < IN_dRightBoundary)
	{
		return true;
	}
	return false;

}

//计算角度的范围
double g_fnAngleRange(double IN_dVar)
{
	double dResult = IN_dVar;
	if (fabs(dResult) > 360)
	{
		int iIN = (int)(IN_dVar);
		dResult = iIN % (360);
	}

	return dResult;
}

//限幅函数
double g_fnIntervalLimit(double IN_dLowBoundary, double IN_dUpBoundary, double IN_dVal)
{
	double dOut = IN_dVal;
	if (IN_dLowBoundary >= IN_dUpBoundary)
	{
		cout << "INPUT BOUNDARY ERROR" << endl;
		return dOut;
	}

	if (dOut > IN_dUpBoundary)
	{
		dOut = IN_dUpBoundary;
	}

	if (dOut < IN_dLowBoundary)
	{
		dOut = IN_dLowBoundary;
	}

	return dOut;
}


void g_fnCalLLAGrid(double IN_dLon_L, double IN_dLon_R, double IN_dLat_D, double IN_dLat_U, C_Vector3 IN_vLLACenter, double IN_dLon_Step, double IN_dLat_Step, int IN_iPointInterval, const vector<C_Vector3>& IN_p_vec_vLLALine, vector<C_Vector3>* Out_p_vec_vLLAGrid)
{
	//未编写
}

//指数函数
double g_fnIntPow(double IN_dVal, int IN_iPow)
{
	double dOut = 1;
	for (int i = 0; i < IN_iPow; i++)
	{
		dOut *= IN_dVal;
	}
	return dOut;
}

//计算高斯噪声
double g_fnGaussNoise(double IN_dMiu, double IN_dSegma, double IN_dr)
{
	int m;
	double s, w, v, ra;
	double r = IN_dr;
	s = 65536.0;
	w = 2053.0;
	v = 13849.0;
	ra = 0.0;
	for (int i = 1; i <= 12; i++)
	{
		r = (r)*w + v;
		m = (int)(r / s);
		r = r - m * s;
		ra = ra + (r) / s;
	}

	ra = IN_dMiu + IN_dSegma * (ra - 6.0);
	return ra;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//插值函数
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//一维内插
double interpolation_inside_1(vector<double> Xs, vector<double> Data, double X)
{
	// Step1：确定Xs向位置
	int id = 0;

	for (int i = 0; i != Xs.size() - 1; i++)
	{
		id = i;
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	double ratio;

	if (X < Xs[0])
	{
		return Data[0];
	}
	else if (X >= Xs[Xs.size() - 1])
	{
		return Data[Xs.size() - 1];
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		return Data[id] + (Data[id + 1] - Data[id]) * ratio;
	}
}

//二维内插
double interpolation_inside_2(vector<double> Xs, vector<double> Ys, vector<double> Data, double X, double Y)
{
	int n_Ys = Ys.size();
	int n_Xs = Xs.size();

	// Step1：确定Xs向位置
	int id = 0;
	for (int i = 0; i != n_Xs - 1; i++)
	{
		id = i;
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	vector<double> Y_Data;
	Y_Data.resize(n_Ys);
	double ratio;

	if (X < Xs[0])
	{
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs];
	}
	else if (X >= Xs[n_Xs - 1])
	{
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + n_Xs - 1];
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + id] + (Data[i * n_Xs + id + 1] - Data[i * n_Xs + id]) * ratio;
	}


	// Step4：Ys向插值
	double	ResultVal;
	if (Y < Ys[0])
	{
		ResultVal = Y_Data[0];
	}

	else if (Y >= Ys[n_Ys - 1])
	{
		ResultVal = Y_Data[n_Ys - 1];
	}
	else
	{
		ratio = (Y - Ys[id]) / (Ys[id + 1] - Ys[id]);
		ResultVal = Y_Data[id] + (Y_Data[id + 1] - Y_Data[id]) * ratio;
	}

	return ResultVal;

}

//三维内插
double interpolation_inside_3(vector<double> Xs, vector<double> Ys, vector<double> Zs, vector<double> Data, double X, double Y, double Z)
{
	int n_Xs = Xs.size();
	int n_Ys = Ys.size();
	int n_Zs = Zs.size();

	// Step1：Xs、Ys向二维插值
	vector<double> Z_Data;
	Z_Data.resize(n_Zs);

	vector<double> data1;
	data1.resize(n_Xs * n_Ys);

	for (int i = 0; i != n_Zs; i++)
	{
		for (int j = 0; j < data1.size(); j++)
		{
			data1[j] = Data[i * n_Xs * n_Ys + j];
		}
		Z_Data[i] = interpolation_inside_2(Xs, Ys, data1, X, Y);
	}

	// Step2：确定Zs向位置
	int id;
	for (int i = 0; i != n_Zs; i++)
	{
		id = i;
		if (Z >= Zs[i] && Z <= Zs[i + 1])
			break;
		else
			continue;
	}

	// Step3：Zs向插值
	double ResultVal;
	double ratio;
	if (Z < Zs[0])
	{
		ResultVal = Z_Data[0];
	}
	else if (Z >= Zs[n_Zs - 1])
	{
		ResultVal = Z_Data[n_Zs - 1];
	}
	else
	{
		ratio = (Z - Zs[id]) / (Zs[id + 1] - Zs[id]);
		ResultVal = Z_Data[id] + (Z_Data[id + 1] - Z_Data[id]) * ratio;
	}

	return ResultVal;
}

//一维外插
double interpolation_outside_1(vector<double> Xs, vector<double> Data, double X)
{
	// Step1：确定Xs向位置
	int id = 0;

	for (int i = 0; i != Xs.size() - 1; i++)
	{
		id = i;
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		return Data[0] - ratio * (Data[1] - Data[0]);
	}
	else if (X >= Xs[Xs.size() - 1])
	{
		ratio = (X - Xs[Xs.size() - 1]) / (Xs[Xs.size() - 1] - Xs[Xs.size() - 2]);
		return Data[Xs.size() - 1] + ratio * (Data[Xs.size() - 1] - Data[Xs.size() - 2]);
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		return Data[id] + (Data[id + 1] - Data[id]) * ratio;
	}
}

//二维外插
double interpolation_outside_2(vector<double> Xs, vector<double> Ys, vector<double> Data, double X, double Y)
{
	int n_Ys = Ys.size();
	int n_Xs = Xs.size();

	// Step1：确定Xs向位置
	int id = 0;
	for (int i = 0; i != n_Xs - 1; i++)
	{
		id = i;
		if (X >= Xs[i] && X <= Xs[i + 1])
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	vector<double> Y_Data;
	Y_Data.resize(n_Ys);
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs] - ratio * (Data[i * n_Xs + 1] - Data[i * n_Xs]);
	}
	else if (X >= Xs[n_Xs - 1])
	{
		ratio = (X - Xs[n_Xs - 1]) / (Xs[n_Xs - 1] - Xs[n_Xs - 2]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + n_Xs - 1] + ratio * (Data[i * n_Xs + n_Xs - 1] - Data[i * n_Xs + n_Xs - 2]);
	}
	else
	{
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + id] + (Data[i * n_Xs + id + 1] - Data[i * n_Xs + id]) * ratio;
	}


	// Step4：Ys向插值
	double	ResultVal;
	if (Y < Ys[0])
	{
		ratio = (Ys[0] - Y) / (Ys[1] - Ys[0]);
		ResultVal = Y_Data[0] - ratio * (Y_Data[1] - Y_Data[0]);
	}

	else if (Y >= Ys[n_Ys - 1])
	{
		ratio = (Y - Ys[n_Ys - 1]) / (Ys[n_Ys - 1] - Ys[n_Ys - 2]);
		ResultVal = Y_Data[n_Ys - 1] + ratio * (Y_Data[n_Ys - 1] - Y_Data[n_Ys - 2]);
	}
	else
	{
		ratio = (Y - Ys[id]) / (Ys[id + 1] - Ys[id]);
		ResultVal = Y_Data[id] + (Y_Data[id + 1] - Y_Data[id]) * ratio;
	}

	return ResultVal;
}

//三维外插
double interpolation_outside_3(vector<double> Xs, vector<double> Ys, vector<double> Zs, vector<double> Data, double X, double Y, double Z)
{
	int n_Xs = Xs.size();
	int n_Ys = Ys.size();
	int n_Zs = Zs.size();

	// Step1：Xs、Ys向二维插值
	vector<double> Z_Data;
	Z_Data.resize(n_Zs);

	vector<double> data1;
	data1.resize(n_Xs * n_Ys);

	for (int i = 0; i != n_Zs; i++)
	{
		for (int j = 0; j < data1.size(); j++)
		{
			data1[j] = Data[i * n_Xs * n_Ys + j];
		}
		Z_Data[i] = interpolation_outside_2(Xs, Ys, data1, X, Y);
	}

	// Step2：确定Zs向位置
	int id;
	for (int i = 0; i != n_Zs; i++)
	{
		id = i;
		if (Z >= Zs[i] && Z <= Zs[i + 1])
			break;
		else
			continue;
	}

	// Step3：Zs向插值
	double ResultVal;
	double ratio;
	if (Z < Zs[0])
	{
		ratio = (Zs[0] - Z) / (Zs[1] - Zs[0]);
		ResultVal = Z_Data[0] - ratio * (Z_Data[1] - Z_Data[0]);
	}
	else if (Z >= Zs[n_Zs - 1])
	{
		ratio = (Z - Zs[n_Zs - 1]) / (Zs[n_Zs - 1] - Zs[n_Zs - 2]);
		ResultVal = Z_Data[n_Zs - 1] + ratio * (Z_Data[n_Zs - 1] - Z_Data[n_Zs - 2]);
	}
	else
	{
		ratio = (Z - Zs[id]) / (Zs[id + 1] - Zs[id]);
		ResultVal = Z_Data[id] + (Z_Data[id + 1] - Z_Data[id]) * ratio;
	}

	return ResultVal;
}

//改进的一维插值外插，插值速度可提升
double interp1_breeze(vector<double> Xs, vector<double> Data, double X)
{
	// Step1：确定Xs向位置
	int id = 0;

	//for (int i = 0; i != n_Xs - 1; i++)
	//{
	//	id = i;
	//	if (X >= Xs[i] && X <= Xs[i + 1])
	//		break;
	//	else
	//		continue;
	//}

	// Step2：Xs向插值
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		return Data[0] - ratio * (Data[1] - Data[0]);
	}
	else if (X >= Xs[Xs.size() - 1])
	{
		ratio = (X - Xs[Xs.size() - 1]) / (Xs[Xs.size() - 1] - Xs[Xs.size() - 2]);
		return Data[Xs.size() - 1] + ratio * (Data[Xs.size() - 1] - Data[Xs.size() - 2]);
	}
	else
	{
		id = upper_bound(Xs.begin(), Xs.end(), X) - Xs.begin() - 1;
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		return Data[id] + (Data[id + 1] - Data[id]) * ratio;
	}

}

//改进的二维插值外插，插值速度可提升
double interp2_breeze(vector<double> Xs, vector<double> Ys, vector<double> Data, double X, double Y)
{
	int n_Ys = Ys.size();
	int n_Xs = Xs.size();

	// Step1：确定Xs向位置
	int id = 0;

	// Step2：Xs向插值
	vector<double> Y_Data;
	Y_Data.resize(n_Ys);
	double ratio;

	if (X < Xs[0])
	{
		ratio = (Xs[0] - X) / (Xs[1] - Xs[0]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs] - ratio * (Data[i * n_Xs + 1] - Data[i * n_Xs]);
	}
	else if (X >= Xs[n_Xs - 1])
	{
		ratio = (X - Xs[n_Xs - 1]) / (Xs[n_Xs - 1] - Xs[n_Xs - 2]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + n_Xs - 1] + ratio * (Data[i * n_Xs + n_Xs - 1] - Data[i * n_Xs + n_Xs - 2]);
	}
	else
	{
		id = upper_bound(Xs.begin(), Xs.end(), X) - Xs.begin() - 1;
		ratio = (X - Xs[id]) / (Xs[id + 1] - Xs[id]);
		for (int i = 0; i != n_Ys; i++)
			Y_Data[i] = Data[i * n_Xs + id] + (Data[i * n_Xs + id + 1] - Data[i * n_Xs + id]) * ratio;
	}


	// Step4：Ys向插值
	double	ResultVal;
	if (Y < Ys[0])
	{
		ratio = (Ys[0] - Y) / (Ys[1] - Ys[0]);
		ResultVal = Y_Data[0] - ratio * (Y_Data[1] - Y_Data[0]);
	}

	else if (Y >= Ys[n_Ys - 1])
	{
		ratio = (Y - Ys[n_Ys - 1]) / (Ys[n_Ys - 1] - Ys[n_Ys - 2]);
		ResultVal = Y_Data[n_Ys - 1] + ratio * (Y_Data[n_Ys - 1] - Y_Data[n_Ys - 2]);
	}
	else
	{
		id = upper_bound(Ys.begin(), Ys.end(), Y) - Ys.begin() - 1;
		ratio = (Y - Ys[id]) / (Ys[id + 1] - Ys[id]);
		ResultVal = Y_Data[id] + (Y_Data[id + 1] - Y_Data[id]) * ratio;
	}

	return ResultVal;
}

//改进的三维插值外插，插值速度可提升
double interp3_breeze(vector<double> Xs, vector<double> Ys, vector<double> Zs, vector<double> Data, double X, double Y, double Z)
{
	int n_Xs = Xs.size();
	int n_Ys = Ys.size();
	int n_Zs = Zs.size();

	// Step1：Xs、Ys向二维插值
	vector<double> Z_Data;
	Z_Data.resize(n_Zs);

	vector<double> data1;
	data1.resize(n_Xs * n_Ys);

	for (int i = 0; i != n_Zs; i++)
	{
		for (int j = 0; j < data1.size(); j++)
		{
			data1[j] = Data[i * n_Xs * n_Ys + j];
		}
		Z_Data[i] = interp2_breeze(Xs, Ys, data1, X, Y);
	}


	// Step3：Zs向插值
	double ResultVal;
	double ratio;
	if (Z < Zs[0])
	{
		ratio = (Zs[0] - Z) / (Zs[1] - Zs[0]);
		ResultVal = Z_Data[0] - ratio * (Z_Data[1] - Z_Data[0]);
	}
	else if (Z >= Zs[n_Zs - 1])
	{
		ratio = (Z - Zs[n_Zs - 1]) / (Zs[n_Zs - 1] - Zs[n_Zs - 2]);
		ResultVal = Z_Data[n_Zs - 1] + ratio * (Z_Data[n_Zs - 1] - Z_Data[n_Zs - 2]);
	}

	else
	{
		int id = upper_bound(Zs.begin(), Zs.end(), Z) - Zs.begin() - 1;
		ratio = (Z - Zs[id]) / (Zs[id + 1] - Zs[id]);
		ResultVal = Z_Data[id] + (Z_Data[id + 1] - Z_Data[id]) * ratio;
	}

	return ResultVal;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*有问题的函数部分*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
//void g_fnCalLLAGrid(double IN_dLon_L, double IN_dLon_R, double IN_dLat_D, double IN_dLat_U, C_Vector3 IN_vLLACenter, double IN_dLon_Step, double IN_dLat_Step, int IN_iPointInterval, const vector<C_Vector3>& IN_p_vec_vLLALine, vector<C_Vector3>* Out_p_vec_vLLAGrid)
	//{
	//	Out_p_vec_vLLAGrid->clear();
	//	if (IN_p_vec_vLLALine.size() == 0)
	//	{
	//		cout << "input traj size zero" << endl;
	//		return;
	//	}
	//	//计算中心点周围的经纬度最大范围
	//	double dMinLon = IN_vLLACenter.x + IN_dLon_L;
	//	double dMaxLon = IN_vLLACenter.x + IN_dLon_R;
	//	double dMinLat = IN_vLLACenter.y + IN_dLat_D;
	//	double dMaxLat = IN_vLLACenter.y + IN_dLat_U;
	//	//判断轨迹是否在范围内
	//	bool bIfExistIn = false;//输入的轨迹是否能与中心点附近范国相交
	//	bool bIfLonAllIn = true; //输入的LLA轨迹经度是否全在中心点周围
	//	bool bIfLatAllIn = true; //输入的LLA轨迹纬度是否全在中心点周围
	//	for (unsigned int i = 0; i < IN_p_vec_vLLALine.size(); i++)
	//	{
	//		if (!bIfExistIn)
	//		{
	//			if ((IN_p_vec_vLLALine[i].x < dMaxLon && IN_p_vec_vLLALine[i].x > dMinLon) && (IN_p_vec_vLLALine[i].y< dMaxLat && IN_p_vec_vLLALine[i].y > dMinLat))
	//			{
	//				bIfExistIn = true;
	//			}
	//			if (bIfLonAllIn)
	//			{
	//				if (IN_p_vec_vLLALine[i].x > dMaxLon || IN_p_vec_vLLALine[i].x < dMinLon)
	//				{
	//					//经度存在超出范围点
	//					bIfLonAllIn = false;
	//				}
	//			}
	//			if (bIfLatAllIn)
	//			{
	//				if (IN_p_vec_vLLALine[i].y > dMaxLat || IN_p_vec_vLLALine[i].y < dMinLat)
	//				{
	//					//经度存在超出范围点
	//					bIfLatAllIn = false;
	//				}
	//			}
	//		}
	//		//分情况判断使用经度还是纬度画框，经度纬度中有一个平行于输入的LLA轨迹，有一个按照正常斗行
	//		int iLonInNum = 0;//经度入点数组号
	//		int iLonOutNum = 0;//经度出点数组号
	//		int iLatInNum = 0;//纬度入点数组号
	//		int iLatOutNum = 0;//纬度出点数组号
	//		if (bIfLatAllIn && !bIfLonAllIn)
	//		{
	//			//经度不都在范国内，纬度都在范国内，按照经度取局部
	//			//搜素经度入点出点
	//			bool bIfStartIn = false;//当前遍历点是否已进入
	//			for (unsigned int i = 0; i < IN_p_vec_vLLALine.size(); i + +)
	//			{
	//				//找入点
	//				if (!bIfStartIn)
	//				{
	//					if (IN_p_vec_vLLALine[i].x > dMinLon && IN_p_vec_vLLALine[i].x < dMaxLon)
	//					{
	//						iLonInNum = i;
	//						bIfStartIn = true;
	//					}
	//				}
	//				//找出点
	//				if (bIfStartIn)
	//				{
	//					if (IN_p_vec_vLLALine[i].x > dMinLon && IN_p_vec_vLLALine[i].x < dMaxLon)
	//					{
	//						iLonOutNum = i;
	//					}
	//					else
	//					{
	//						break;
	//					}
	//				}
	//			}

	//			/// 按照输入的间隔取经度，按照输入的纬度步长和上下限取纬度，生成网格
	//			for (int i = iLonInNum; i <= iLonOutNum; i++)
	//			{
	//				if (i % IN_iPointInterval != 0 && i <= iLonOutNum)
	//				{
	//					continue;
	//				}
	//				double dLoopLat = IN_p_vec_vLLALine[i].y + IN_dLat_D;
	//				while (dLoopLat <= IN_p_vec_vLLALine[i].y + IN_dLat_U + SMALLNUM)
	//				{
	//					C_Vector3 vLLAGrid;
	//					vLLAGrid.SetElems(IN_p_vec_vLLALine[i].x, dLoopLat, 0);
	//					Out_p_vec_vLLAGrid->push_back(vLLAGrid);

	//					dLoopLat += IN_dLat_Step;
	//				}
	//			}
	//		}
	//		else if (bIfLonAllIn && !bIfLatAllIn)
	//		{
	//			//经度都在范围内，纬度不都在
	//			//搜索纬度入点出点
	//			bool bIfStartIn = false;//当前遍历点是否已进入
	//			for (unsigned int i = 0; i < IN_p_vec_vLLALine.size(); i++)
	//			{
	//				//找入点
	//				if (!bIfStartIn)
	//				{
	//					if (IN_p_vec_vLLALine[i].y > dMinLat && IN_p_vec_vLLALine[i].y < dMaxLat)
	//					{
	//						iLatInNum = i;
	//						bIfStartIn = true;
	//					}
	//				}
	//				// 找出点
	//				if (bIfStartIn)
	//				{
	//					if (IN_p_vec_vLLALine[i].y > dMinLat && IN_p_vec_vLLALine[i].y < dMaxLat)
	//					{
	//						iLatOutNum = i;
	//					}
	//					else
	//					{
	//						break;
	//					}
	//				}
	//			}

	//			// 按照输入的间隔取纬度，按照输入的经度步长和上下限取经度，生成网格
	//			for (int i = iLatInNum; i<= iLatOutNum; i++)
	//			{
	//				if (i % IN_iPointInterval != 0 && i != iLatOutNum)
	//				{
	//					continue;
	//				}
	//				double dLoopLon = IN_p_vec_vLLALine[i].x +IN_dLon_L;
	//				while (dLoopLon <= IN_p_vec_vLLALine[i].x +IN_dLon_R + SMALLNUM)
	//				{
	//					C_Vector3 vLLAGrid;
	//					vLLAGrid.SetElems(dLoopLon, IN_p_vec_vLLALine[i].y, 0);
	//					Out_p_vec_vLLAGrid ->push_back(vLLAGrid);
	//					dLoopLon += IN_dLon_Step;
	//				}
	//				else if (bIfLonAllIn && bIfLatAllIn)
	//				{
	//					//所有经纬度都在范围内
	//					//看经纬度哪个变化的角度更大（防止出现乖直的情况）
	//					double dLonChange = fabs(IN_p_vec_vLLALine[0].x - IN_p_vec_vLLALine.back().x);
	//					double dLatChange = fabs(IN_p_vec_vLLALine[0].y - IN_p_vec_vLLALine.back().y);
	//					//比较
	//					if (dLonChange >= dLatChange)
	//					{

	//						//经度变化大，纬度可能出现不变情况，用经度框纬度
	//						for(unsigned int i = 0; i < IN_p_vec_vLLALine.size(); i++)
	//						{
	//							if (i % IN_iPointInterval != 0 && i != IN_p_vec_vLLALine.size() - 1)
	//							{
	//								continue;
	//							}
	//							double dLoopLat = IN_p_vec_vLLALine[i].y + IN_dLat_D;
	//							while (dLoopLat <= IN_p_vec_vLLALine[i].y + IN_dLat_U + SMALLNUM)
	//							{
	//								C_Vector3 vLLAGrid;
	//								vLLAGrid.SetElems(IN_p_vec_vLLALine[i].x, dLoopLat, 0);
	//								Out_p_vec_vLLAGrid ->push_back(vLLAGrid);
	//								dLoopLat += IN_dLat_Step;
	//							if (dLonChange < dLatChange)
	//							{

	//								// 纬度变化大，经度可能出现不变情况。用纬度框经度
	//					for (unsigned int i = 0 : i < INp_vec_VLLALine.size0; i + t)

	//		
	//}



	//void g_fnTranTarEarthToLau(C_Vector3 IN_vLauLLA, double IN_dLauAzi, const C_Vector3 & IN_vEarthTarPose, const C_Vector3 & IN_vEarthTarVelo, C_Vector3 * OUT_vLauTarPose, C_Vector3 * OUT_vLauTarVelo)
	//{
	//	//计算初始地心纬度角、初始地理地心纬度角差值
	//	double dEarthPhi0 = TargetPred::B2phi(IN_vLauLLA.y * DEG2RAD);
	//	double dLauMu0 = IN_vLauLLA.y * DEG2RAD - dEarthPhi0;
	//	double dLauR0 = EARTH_RAVG + IN_vLauLLA.z;
	//	// 计算发射点地心矢径在发射系下矢量
	//	C_Vector3 vLauR0;
	//	vLauR0.x = -dLauR0 * sin(dLauMu0) * cos(IN_dLauAzi * DEG2RAD);
	//	vLauR0.y = dLauR0 * cos(dLauMu0);
	//	vLauR0.z = dLauR0 * sin(dLauMu0) * sin(IN_dLauAzi * DEG2RAD);
	//	// 计算目标地心矢径在发射系下分量
	//	C_Vector3 vLauTarEarthPose;
	//	vLauTarEarthPose = getMatrix_EarthToLaunch(IN_vLauLLA.x * DEG2RAD, IN_vLauLLA.y * DEG2RAD, IN_dLauAzi * DEG2RAD) * IN_vEarthTarPose;
	//	//计算目标在发射系下位置
	//	*OUT_vLauTarPose = vLauTarEarthPose - vLauR0;
	//	// 计算目标在发射系下位置
	//	*OUT_vLauTarVelo = getMatrix_EarthToLaunch(IN_vLauLLA.x * DEG2RAD, IN_vLauLLA.y * DEG2RAD, IN_dLauAzi * DEG2RAD) * IN_vEarthTarVelo;
	//	
	//	return;
	//}
	*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////













