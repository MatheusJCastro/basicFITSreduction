<style>
body {
  text-align: justify;
}

h1, h2, h3, h4, h5, h6 {
  font-weight: bold;
}

img {
  display: block;
  margin-left: auto;
  margin-right: auto;
}

.center {
  display: flex;
  justify-content: center;
  align-items: center;
}

code, .language-bash {color: #cc992c;
font-weight: bold;}
</style>

# How to Reduce FITS Images - Basic Steps

*Made by: Matheus J. Castro*  
*Email: mcastro@lsw.uni-heidelberg.de*  

*GitHub Copilot was used to assist in writing this document.*

## Table of Contents

<hr/>

- [How to Reduce FITS Images - Basic Steps](#how-to-reduce-fits-images---basic-steps)
  - [Table of Contents](#table-of-contents)
  - [Bias Reduction](#bias-reduction)
  - [Flat Reduction](#flat-reduction)
  - [The final Equation](#the-final-equation)
  - [RGB Images](#rgb-images)
  - [The Code](#the-code)

## Bias Reduction

<hr/>

The bias is the electronic offset of the CCD camera. It is a small signal that is added to every pixel in the image to prevent negative values. To reduce the bias, you need to take several bias frames (images taken with zero exposure time) and combine them to create a master bias frame. This master bias frame can then be subtracted from your science images.

The bias frames can be taken at any time, as long as the camera is in the same state as when the science images were taken. It is recommended to take at least 10 bias frames to create a good master bias frame.

The equation for bias reduction is:

```math
Master\_Bias = median(Bias\_Frames) \\
Reduced\_Image = Science\_Image - Master\_Bias
```

## Flat Reduction

<hr/>

The flat is the response of the CCD camera to a uniform light source. It is used to correct for variations in the sensitivity of the camera. To reduce the flat, you need to take several flat frames (images taken with a uniform light source) and combine them to create a master flat frame. This master flat frame can then be used to normalize your science images.

The flat frames should be taken with the same filter and focus as the science images. It is recommended to take at least 10 flat frames to create a good master flat frame.

Sky flats can be taken during twilight, while dome flats can be taken using a light source inside the observatory dome. It is important to ensure that the flat frames are not saturated and that they cover the entire field of view of the camera.

Flats should also have the bias subtracted before combining them to create the master flat frame. The equation for flat reduction is:

```math
Master\_Flat = median(Flat\_Frames - Master\_Bias)
```

The flat should also be normalized to have a median value of 1. This can be done by dividing the master flat by its median value:

```math
Normalized\_Master\_Flat = \frac{Master\_Flat}{median(Master\_Flat)}
```

## The final Equation

<hr/>

The equation which includes both flat and bias reduction is:

```math
Reduced\_Image = \frac{Science\_Image - Master\_Bias}{ \frac{Master\_Flat}{median(Master\_Flat)}}
```

## RGB Images

<hr/>

Astronomical images are often taken in different filters, such as red, green, and blue. Each filter captures a different wavelength of light, but the images are still black and white. By combining them you can create a colored image.  

To create a color image, you need to reduce the images taken in each filter separately and then combine them into a single RGB image. It is not necessary to use the R, G, and B filters specifically; you can use any combination of filters to create a color image, and in different orders.

Different filter will show different features of the object being imaged. For example, using a narrowband filter can highlight specific emission lines, while using a broadband filter can show the overall structure of the object. 

- An H-alpha filter is a narrowband filter that isolates hydrogen's H-alpha emission line at 6563 Angstroms, making it useful for imaging emission nebulae, including those associated with the ionized gas in star-forming regions. 
- OIII, at ~5000 Angstroms, and SII, at ~6720 Angstroms, are narrowband filters that isolate emission lines from doubly ionized oxygen and singly ionized sulfur, respectively. They are commonly used to reveal structures in planetary nebulae and supernova remnants..

Before combining the images, it is important to ensure that they are aligned.

## The Code

<hr/>

A simple code to reduce FITS images is available in [here](https://github.com/MatheusJCastro/basicFITSreduction). You can install it using pip:

```bash
pip install git+https://github.com/MatheusJCastro/basicFITSreduction.git
```

To use the code, the data folder structure should be as follows:

```bash
Root_Folder/
├── Data_Folder/
│   ├── ObjectName-001-R.fits
│   ├── ObjectName-002-R.fits
│   ├── ObjectName-001-G.fits
│   ├── ObjectName-002-G.fits
│   └── ...
├── bias-001.fits
├── bias-002.fits
├── flat-001-R.fits
├── flat-002-R.fits
├── flat-001-G.fits
└── ...
```

A folder named `calibrated` will be created inside the `Data_Folder` to store the reduced images. 

To run the code, check the [Example.py](https://github.com/MatheusJCastro/basicFITSreduction/blob/main/Example.py) in the GitHub repository.

Obs: for a successful reduction, the name of the filters should be written in the same way in every file. 