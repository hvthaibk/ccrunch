# ccrunch

Welcome to gEMtools – a collection of programs for processing cryo-electron microscopy ("cryo-EM") datasets. The overall aim in cryo-EM is to reconstruct the 3D structures of large macromolecules from multiple 2D images, or "micrographs", which have been captured by an electron microscope. Because many of the steps involved in processing such 2D micrographs and the resulting 3D volumes are computationally expensive, the gEMtools suite uses fast Fourier transform (FFT) correlation techniques to accelerate the calculations, and it uses parallel processing techniques to distribute the FFT calculations over multiple CPU cores in a workstation and multiple nodes in a computer cluster. A particular novelty of the gEMtools programs is that they can use multiple graphics processor units (GPUs) to accelerate the correlation calculations even further.

Related information: http://gem.loria.fr/
