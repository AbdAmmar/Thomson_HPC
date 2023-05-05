
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
      <a href="#prerequisite">Getting Started</a>
      <ul>
        <li><a href="#prerequisite">Prerequisites</a></li>
        <li><a href="#build">Build</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#example">Example</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

# Thomson
Thomson is program to solve thomson problem in torus in one,two and three dimension based on FORTAN90

![Screenshot](src/GUI.png)

# Prerequisite

* 1 - Fortran compiler (gfortran recommended):

```
sudo apt install gfortran
```
* 2 - Fortran library (llapack, openmp):
  
```
sudo apt-get install liblapack-dev libopenblas-dev libomp-dev
```
* 3 - For animation (gnuplot, FFmpeg):

```
sudo apt-get install gnuplot ffmpeg
```
* 4 - Python (Python 3, not working using python 2 because it need some library )

```
sudo apt-get install python3 
```
* 5 - Python library (numpy,tkinter,customtkinter) to use the python GUI

```
pip3 install customtkinter tkinter numpy 
```

 
# Build

* build the program
  
```sh
make Thomson
```
* clean the bin files
```sh
make clean
```
# Usage 

* Using Python code

```sh
python Thomson.py
``` 
* Thomson

```sh
./Thomson "The input file name"
```
  
# Example

there are some example of the inputfile in the [examples folder]() you can try to use them.
  
## Contact

Stefano Evangelisti - stefano.lcpq@gmail.com

Amer Alrakik - alrakikamer@gmail.com

Project Link: [https://github.com/ALRAKIK/Thomson](https://github.com/ALRAKIK/Thomson)

## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

## Acknowledgments

Special thanks go to Prof.Arjan Berger and Rony J. Letona for your ideas and help.
