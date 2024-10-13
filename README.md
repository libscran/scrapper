# R bindings to C++ single-cell analysis libraries

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/scrapper.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/scrapper.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/scrapper/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/scrapper.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/scrapper.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/scrapper/)|

The **scrapper** package implements R bindings to C++ code for analyzing single-cell (expression) data, mostly various [**libscran**](https://github.com/libscran) libraries. 
Each function performs an individual step in the single-cell analysis workflow, ranging from quality control to clustering and marker detection.
It is mostly intended for other Bioconductor package developers to build more user-friendly end-to-end workflows.
