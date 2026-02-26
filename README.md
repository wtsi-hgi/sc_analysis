<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">Single Cell Multiome Analysis</h3>

  <p align="center">
    TF-induced expresssion profile
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
Synthetic Lineage Project -- by Connor Rogerson

To learn sequence to expression rules, we need to figure out how TFs interact with DNA and affect gene expression.
<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started
The pipeline is designed for single-cell multiome analysis, including scRNA-seq, scATAC-seq, and transcription factor profiling.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- PREREQUISITES -->
### Prerequisites
<details>
<summary><b>Python Libraries</b></summary>

    pandas = 2.2.2
    numpy = 2.1.2
    biopython = 1.85
    polars = 1.25.2
</details>

<details>
<summary><b>R Packages</b></summary>

    optparse = 1.7.4
    tidyverse = 2.0.0
    data.table = 1.15.4
    ggplot2 = 3.5.2
    ggrepel = 0.9.5
    corrplot = 0.92
    GenomicRanges = 1.54.1
    Seurat = 5.3.0
    Signac = 1.12.0
    SeuratWrappers = 0.3.4
    SoupX = 1.6.2
    DoubletFinder = 2.0.6
    GenomeInfoDb = 1.38.8
</details>

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

### Sample Sheet


<p align="right">(<a href="#top">back to top</a>)</p>

<!-- LICENSE -->
## License
Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTACT -->
## Contact
Example

<p align="right">(<a href="#top">back to top</a>)</p>

TF plot:
y -> median no. of TFs
x -> differnt tf count filters (> 1,2,3 etc..)

report:
summary/stats -> html report

h5 outputs + matrix tvs files

