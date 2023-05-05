
<br />
<div align="center">
  <a href="https://https://github.com/ALRAKIK/Thomson">
    <img src="src/logo.png" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">Thomson Project</h3>

  <p align="center">
    Solve Thomson problem in tours in one, two and three dimension
    <br />
    <a href="https://github.com/ALRAKIK/Thomson"><strong>Explore the Project »</strong></a>
    <br />
    <br />
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#thomson">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisite">Prerequisites</a></li>
        <li><a href="#build">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

# Thomson
Thomson is program to solve thomson problem in torus in one,two and three dimension based on FORTAN90

![Screenshot](src/GUI.png)

# Prerequisite
* 1 - Fortran compiler (gfortran recommended):

  ```sh
  sudo apt install gfortran
  ```
* 2 - Fortran library (llapack, openmp):
  
  ```sh
  sudo apt-get install liblapack-dev libopenblas-dev libomp-dev
  ```
* 3 - For animation (gnuplot, FFmpeg):

  ```sh
  sudo apt-get install gnuplot ffmpeg
  ```
* 4 - Python (Python 3, not working using python 2 because it need some library )

 ```sh
  sudo apt-get install python3 
 ```
 * 5 - Python library using pip

 ```sh
  sudo apt-get install python3 
 ```







  
# Build

* build the program
  
  ```sh
  make Thomson
  ```
 


# Run 

* Using Python code

  ```sh
   python Thomson.py
  ``` 
* Thomson

  ```sh
  ./Thomson "The input file name"
  ```
  
* there are some example of the inputfile in the example folder you can use anyone of them or make your input
  
  
<!-- CONTACT -->
## Contact

Stefano Evangelisti - stefano.lcpq@gmail.com

Amer Alrakik - alrakikamer@gmail.com

Project Link: [https://github.com/ALRAKIK/Thomson](https://github.com/ALRAKIK/Thomson)

## License

Distributed under the MIT License. See `LICENSE.txt` for more information.
