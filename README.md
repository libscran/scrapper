# R bindings to C++ single-cell analysis libraries

The **scrapper** package implements R bindings to C++ code for analyzing single-cell (expression) data, mostly various [**libscran**](https://github.com/libscran) libraries. 
Each function performs an individual step in the single-cell analysis workflow, ranging from quality control to clustering and marker detection.
It is mostly intended for other Bioconductor package developers to build more user-friendly end-to-end workflows.
