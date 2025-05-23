#include <Python.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <iostream>
//#include <numpy/arrayobject.h>
//#include <pybind11/pybind11.h>
using std::vector;
using std::logic_error;
using std::string;
//namespace py = pybind11;

void printS(string input, char endl = '\n') {
    std::cout << input << endl;
}

static float mean(vector<float> data) {
    float sum = 0;
    for (int i = 0; i < data.size(); i++) {
        sum += data[i];
    }
    return sum / data.size();
}

static vector<float> allPos(vector<float> data) {
    vector<float> posdata;
    for (int i = 0; i < data.size(); i++) {
        if (data[i] >= 0) {
            posdata.push_back(data[i]);
        }
    }
    return posdata;
}

static float median(vector<float> data) {
    float medianValue;
    std::sort(data.begin(), data.end());

    if (data.size() % 2 == 0) {
        medianValue = (data[data.size() / 2] + (data[(data.size() / 2) - 1])) / 2;
    }
    else {
        medianValue = data[data.size() / 2];
    }
    return medianValue;
}

static double stdev(vector<float> data) {
    double sum = 0;
    float meanValue = mean(data);
    for (int i = 0; i < data.size(); i++) {
        sum += std::pow(data[i] - meanValue, 2);
    }
    sum = sum / data.size();
    return std::pow(sum, 0.5);
}

static vector<float> listToVector_Float(PyObject* incoming) {
    vector<float> data;
    if (PyList_Check(incoming)) {
        for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
            PyObject* value = PyList_GetItem(incoming, i);
            data.push_back(PyFloat_AsDouble(value));
        }
    }
    else {
        throw logic_error("Passed PyObject pointer was not a list!");
    }
    return data;
}

static vector<int> listToVector_Int(PyObject* incoming) {
    vector<int> data;
    if (PyList_Check(incoming)) {
        for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
            PyObject* value = PyList_GetItem(incoming, i);
            data.push_back(PyLong_AsLong(value));
        }
    }
    else {
        throw logic_error("Passed PyObject pointer was not a list!");
    }
    
    return data;
}

static vector<vector<float>> listToVector_Float2d(PyObject* incoming) {
    vector<vector<float>> data;
    if (PyList_Check(incoming)) {
        for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
            PyObject* value = PyList_GetItem(incoming, i);
            vector<float> vecvalue = listToVector_Float(value);
            data.push_back(vecvalue);
        }
    }
    else {
        throw logic_error("Passed PyObject pointer was not a list or tuple!");
    }
    return data;
}

vector<vector<int>> listToVector_int2d(PyObject* incoming) {
    vector<vector<int>> data;
    if (PyList_Check(incoming)) {
        for (Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
            PyObject* value = PyList_GetItem(incoming, i);
            vector<int> vecvalue = listToVector_Int(value);
            data.push_back(vecvalue);
        }
    }
    else {
        throw logic_error("Passed PyObject pointer was not a list or tuple!");
    }
    return data;
}

PyObject* vectorToList_Float2d(const vector<vector<float>>& data) {
    PyObject* listObj = PyList_New(data.size());
    if (!listObj) throw logic_error("Unable to allocate memory for Python list");
    for (unsigned int i = 0; i < data.size(); i++) {
        PyObject* templistObj = PyList_New(data[i].size());
        for (unsigned int j = 0; j < data[i].size(); j++) {
            PyObject* num = PyFloat_FromDouble((double)data[i][j]);
            if (!num) {
                Py_DECREF(listObj);
                throw logic_error("Unable to allocate memory for Python list");
            }
            PyList_SET_ITEM(templistObj, j, num);
        }
        PyList_SET_ITEM(listObj, i, templistObj);
    }
    return listObj;
}

PyObject* vectorToList_int2d(const vector<vector<int>>& data) {
    PyObject* listObj = PyList_New(data.size());
    if (!listObj) throw logic_error("Unable to allocate memory for Python list");
    for (unsigned int i = 0; i < data.size(); i++) {
        PyObject* templistObj = PyList_New(data[i].size());
        for (unsigned int j = 0; j < data[i].size(); j++) {
            PyObject* num = PyLong_FromLong((int)data[i][j]);
            if (!num) {
                Py_DECREF(listObj);
                throw logic_error("Unable to allocate memory for Python list");
            }
            PyList_SET_ITEM(templistObj, j, num);
        }
        PyList_SET_ITEM(listObj, i, templistObj);
    }
    return listObj;
}

vector<vector<int>> generateMask(vector<vector<float>> dataArray, vector<vector<int>> basemask,  vector<vector<int>> binArray,
        uint16_t nbins,  float stdevs, uint16_t threshold) {
    vector<vector<float>> binnedData(nbins);
    vector<vector<vector<int>>> indexes(nbins);

    if ( ! ((dataArray.size() == basemask.size() && dataArray.size() == binArray.size()) &&
        (dataArray[0].size() == basemask[0].size() && dataArray[0].size() == binArray[0].size()))) {
        std::cout << "array size mismatch\n";
        std::cout << "data shape: " << dataArray.size() << ',' << dataArray[0].size() << '\n';
        std::cout << "mask shape: " << basemask.size() << ',' << basemask[0].size() << '\n';
        std::cout << "bin shape: " << binArray.size() << ',' << binArray[0].size() << '\n';
        PyErr_BadArgument();
        //throw std::runtime_error("array mismatch");
    }
    int y;
    int x;
    for ( y = 0; y < dataArray.size(); y++) {
        for ( x = 0; x < dataArray[0].size(); x++) {
            if ((dataArray[y][x] < 0) && (basemask[y][x]) == 0){
                basemask[y][x] = 1;
                continue;
            }
            if ((dataArray[y][x] >= 0) && (basemask[y][x] == 0)) {
                binnedData[binArray[y][x]].push_back(dataArray[y][x]);
                indexes[binArray[y][x]].push_back({ y,x });
            }
        }
    }
    float mediandata;
    double stdevData;
    for (int i = 0; i < nbins; i++) {
        vector<float> data = binnedData[i];
        if (data.size() == 0) {
            continue;
        }
        mediandata = median(data);
        stdevData = stdev(data);
        for (int j = 0; j < indexes[i].size(); j++) {
            y = indexes[i][j][0];
            x = indexes[i][j][1];
            if (dataArray[y][x] > mediandata + (stdevData * stdevs) || dataArray[y][x] > mediandata + stdevData + threshold) {
                basemask[y][x] = 1;
            }
        }
    }
    return basemask;
}

PyObject* makeMaskCP(PyObject*, PyObject* argTup) {
    
    if (!PyTuple_Check(argTup)) {
        string message = "Passed PyObject pointer was not a tuple!";
        std::cout << message << '\n';
        throw logic_error(message);
    }
    
    PyObject* dataPy = PyTuple_GetItem(argTup, 0);
    PyObject* baseMaskPy = PyTuple_GetItem(argTup, 1);
    PyObject* binArrayPy = PyTuple_GetItem(argTup, 2);
    PyObject* nbinsPy = PyTuple_GetItem(argTup, 3);
    PyObject* stdevsPy = PyTuple_GetItem(argTup, 4);
    PyObject* thresholdPy = PyTuple_GetItem(argTup, 5);
    

    vector<vector<float>> dataArray = listToVector_Float2d(dataPy);
    vector<vector<int>> baseMask = listToVector_int2d(baseMaskPy);
    vector<vector<int>> binArray = listToVector_int2d(binArrayPy);
    int nbins = PyLong_AsLong(nbinsPy);
    float stdevs = PyFloat_AsDouble(stdevsPy); 
    int threshold = PyLong_AsLong(thresholdPy);
    

    vector<vector<int>> newmask = generateMask(dataArray, baseMask, binArray, nbins, stdevs,
        threshold);
    
    PyObject* newmaskPy = vectorToList_int2d(newmask);
    printS("mask calculated");
    return newmaskPy;
}



static PyMethodDef azimuthCP_methods[] = {
    // The first property is the name exposed to Python, fast_tanh
    // The second is the C++ function with the implementation
    // METH_O means it takes a single PyObject argument
    {"makeMaskCP",(PyCFunction)makeMaskCP,METH_O,nullptr},
    // Terminate the array with an object containing nulls
    { nullptr, nullptr, 0, nullptr }
};


static PyModuleDef azimuthCP_module = {
    PyModuleDef_HEAD_INIT,
    "azimuthCP",                        // Module name to use with Python import statements
    "adds a faster mask generator function",  // Module description
    0,
    azimuthCP_methods                   // Structure that defines the methods of the module
};

PyMODINIT_FUNC PyInit_azimuthCP() {
    return PyModule_Create(&azimuthCP_module);
}