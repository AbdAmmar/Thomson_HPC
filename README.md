
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

# Thomson
Thomson is program to solve thomson problem in torus in one,two and three dimension based on FORTAN90

![Screenshot](src/GUI.png)

# Prerequisite
* compiler (gfortran recommended):

  ```sh
  sudo apt install gfortran
  ```
* library (llapack, openmp):
  
  ```sh
  sudo apt-get install liblapack-dev libopenblas-dev libomp-dev
  ```
* for animation (gnuplot, FFmpeg):

  ```sh
  sudo apt-get install gnuplot ffmpeg
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
