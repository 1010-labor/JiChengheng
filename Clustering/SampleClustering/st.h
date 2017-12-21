#pragma once
// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include <stdio.h>
#include <tchar.h>

#include "math.h"
#include <list>
#include "iomanip"
#include "set"
#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include <iterator>
#include <fstream>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <tchar.h>
#include "assert.h"

using namespace std;

class Record
{
public:
	Record() {};
	Record(const Record& copy,int i) {
		this->index = copy.index;
		this->_x = copy._x;
		this->_y = copy._y;
		this->index = i;
	}
	int index;
	int isBorder;
	double _x;
	double _y;
	
};
struct dist {
	int with;
	double dis;
};
struct dens {
	int index;
	double density;
};
