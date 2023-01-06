#include <python3.8/Python.h>
#include <bits/stdc++.h>
#include "Param.h"
#include "CorrectDetector.h"
#include "Common/Mmap.h"
#include "XSketch.h"
using namespace std;

extern "C"
{
    XSketch<uint64_t>* XSketch_new(uint32_t memory, double var_input, double error_input, double ratio_input, int CellNum, int S1Wins, double potential) {
        init_matrix();
        return new XSketch<uint64_t>(memory, var_input, error_input, ratio_input, CellNum, S1Wins, potential);
    }

    void insert(XSketch<uint64_t>* sm, uint64_t id, uint32_t timestamp) {
        sm->insert(id, timestamp);
        return;
    }
    
    PyObject* query(XSketch<uint64_t>* sm) {
        std::vector<std::pair<uint64_t, uint32_t> > ret = sm->query();
        // std::cout << "Query size: " << ret.size() << std::endl;
        PyObject* listObj = PyList_New(ret.size() * 2);
        for (uint32_t i = 0; i < ret.size(); ++i) {
            PyObject *x = PyLong_FromLong((uint64_t)ret[i].first);
            PyObject *y = PyLong_FromLong((uint32_t)ret[i].second);
            PyList_SET_ITEM(listObj, i * 2, x);
            PyList_SET_ITEM(listObj, i * 2 + 1, y);
        }
        return listObj;
    }

    PyObject* report(XSketch<uint64_t>* sm) {
        std::vector<Report_Slot<uint64_t> > ret = sm->report();
        // std::cout << "Report size: " << ret.size() << std::endl;
        PyObject* listObj = PyList_New(ret.size() * 3);
        for (uint32_t i = 0; i < ret.size(); ++i) {
            PyObject *x = PyLong_FromLong((uint64_t)ret[i].id);
            PyObject *y = PyLong_FromLong((uint32_t)ret[i].start_window);
            PyObject *z = PyLong_FromLong((uint32_t)ret[i].end_window);
            PyList_SET_ITEM(listObj, i * 3, x);
            PyList_SET_ITEM(listObj, i * 3 + 1, y);
            PyList_SET_ITEM(listObj, i * 3 + 2, z);
        }
        return listObj; 
    }

    void eval(XSketch<uint64_t>* sm, double var_input, double error_input, uint32_t run_length) {

        std::string PATH = "/share/datasets/CAIDA2016/1.dat";
	    LoadResult load_result = Load(PATH.c_str());
        CAIDA_Tuple *dataset = (CAIDA_Tuple*)load_result.start;
        uint32_t length = load_result.length / sizeof(CAIDA_Tuple);

        CorrectDetector<uint64_t>* correct_detector = new CorrectDetector<uint64_t>(error_input, var_input);
	    for (uint32_t i = 0; i < run_length; ++i) {
		    correct_detector->insert(dataset[i%length].id, i);
	    }

        std::map<uint32_t, std::map<uint64_t, uint32_t> > predict = sm->predict();
        std::map<uint64_t, std::map<uint32_t, uint32_t> > history = correct_detector->get_history();

        int correct = 0, total = 0;
        for (auto i : predict) {
            for (auto j : i.second) {
                int predict_result = j.second, truth = history[j.first][i.first];
                if (abs(truth - predict_result) <= 5 || (2 * truth >= predict_result && truth <= 2 * predict_result)) {
                    correct++;
                }
                total++;
            }
        }
        std::cout << "Predict: " << total << " Correct: " << correct << "\n";
        std::cout << "Accuracy: " << std::fixed << std::setprecision(4) << 100.0 * correct / total << "\n";

        return;
    }

}

