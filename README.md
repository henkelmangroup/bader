# Bader Charge Analysis

## News

**08/19/23 - [Version 1.05
Released](http://theory.cm.utexas.edu/forum/viewtopic.php?f=1&t=2&p=9465#p9465)**\
Proper reading of negative direct coordinates. Thanks to Yuri Mastrikov
for identifying this issue.

## Introduction

[Richard Bader](http://www.chemistry.mcmaster.ca/bader/), from McMaster
University, developed an intuitive way of dividing molecules into atoms.
His definition of an atom is based purely on the electronic charge
density. Bader uses what are called zero flux surfaces to divide atoms.
A zero flux surface is a 2-D surface on which the charge density is a
minimum perpendicular to the surface. Typically in molecular systems,
the charge density reaches a minimum between atoms and this is a natural
place to separate atoms from each other.

Besides being an intuitive scheme for visualizing atoms in molecules,
Bader\'s definition is often useful for charge analysis. For example,
the charge enclosed within the Bader volume is a good approximation to
the total electronic charge of an atom. The charge distribution can be
used to determine multipole moments of interacting atoms or molecules.
Bader\'s analysis has also been used to define the hardness of atoms,
which can be used to quantify the cost of removing charge from an atom.

## Program Overview

We have developed a fast algorithm for doing Bader\'s analysis on a
charge density grid. The program (see below) can read in charge
densities in the [VASP](https://www.vasp.at/) CHGCAR format, or the
[Gaussian](https://www.gaussian.com/) CUBE format. The program outputs
the total charge associated with each atom, and the zero flux surfaces
defining the Bader volumes.

## Download

Select the appropriate platform to download a binary of the Bader
analysis program:

- [Linux x86-64](https://github.com/henkelmangroup/bader/releases/latest/download/bader_lnx_64.tar.gz) (ifort)
- [Mac OS X](https://github.com/henkelmangroup/bader/releases/latest/download/bader_osx_gfortran.tar.gz) (gfortran)

The F90 source code is also available:

- [Source Code](https://github.com/henkelmangroup/bader/releases/latest/download/bader.tar.gz) (v1.05 08/19/23)

## Running the Program

The program can be run with the command


        bader chargefile

It will automatically determine if the chargefile is a VASP CHGCAR file
or a Gaussian CUBE file. The only required input argument is the name of
the charge density file.

## Command line arguments and output files

The following options can be used when running the Bader analysis
program.


        bader [ -c bader | voronoi ]
              [ -n bader | voronoi ]
              [ -b neargrid | ongrid | weight ]
              [ -r refine_edge_method ]
              [ -ref reference_charge ]
              [ -vac off | auto | vacuum_density ]
              [ -p all_atom | all_bader ]
              [ -p sel_atom | sel_bader ] [volume list or range ]
              [ -p sum_atom | sum_bader ] [ volume list or range ]
              [ -p atom_index | bader_index ]
              [ -i cube | chgcar ]
              [ -cp ]
              [ -h ] [ -v ]
              chargefile

To get a description of the options, run `bader -h`.

## Output files

The following output files are generated: ACF.dat, BCF.dat,
AtomVolumes.dat.

**ACF.dat** contains the coordinates of each atom, the charge associated
with it according to Bader partitioning, percentage of the whole
according to Bader partitioning and the minimum distance to the surface.
This distance should be compared to maximum cut-off radius for the core
region if pseudo potentials have been used.

**BCF.dat** contains the coordinates of each Bader maxima, the charge
within that volume, the nearest atom and the distance to that atom.

**AtomVolumes.dat** contains the number of each volume that has been
assigned to each atom. These numbers correspond to the number of the
BvAtxxxx.dat files.

The Bader volumes can be written using the print options.


        bader [ -p all_atom | all_bader ] chargefile
        bader [ -p sel_atom | sel_bader ] [ volume list or range ] chargefile
        bader [ -p sum_atom | sum_bader ] [ volume list or range ] chargefile
        bader [ -p atom_index | bader_index ] chargefile

**-p none** The default is to write no charge density files.

**-p all_atom** Combine all volumes associated with an atom and write to
file. This is done for all atoms and written to files named
BvAtxxxx.dat. The volumes associated with atoms are those for which the
maximum in charge density within the volume is closest to the atom.

**-p all_bader** Write all Bader volumes (containing charge above
threshold of 0.0001) to a file. The charge distribution in each volume
is written to a separate file, named Bvolxxxx.dat. It will either be of
a CHGCAR format or a CUBE file format, depending on the format of the
initial charge density file. These files can be quite large, so this
option should be used with caution.

**-p sel_atom** Write the selected atomic volumes, read from the
subsequent list or range of volumes.

**-p sel_bader** Write the selected Bader volumes, read from the
subsequent list or range of volumes.

**-p sum_atom** Write the sum of selected atomic volumes, read from the
subsequent list of volumes.

**-p sum_bader** Write the sum of selected Bader volumes, read from the
subsequent list of volumes.

**-p atom_index** Write the atomic volume index to a charge density
file.

**-p bader_index** Write the Bader volume index to a charge density
file.

## Visualization

The Bader volumes can be written and visualized with the [VASP Data
Viewer](https://vaspview.sourceforge.net/),
[VMD](https://www.ks.uiuc.edu/Research/vmd/),
[Jmol](https://jmol.sourceforge.net),
[VESTA](https://jp-minerals.org/vesta/en/), or a cube file viewer (such
as GaussView) for Gaussian cube files.

## Examples

Examples can be found in the \'examples\' directory or downloaded here:


- [NaCl](https://raw.githubusercontent.com/henkelmangroup/bader/main/examples/NaCl.tar.gz) crystal (vasp chgcar)
- [NaCl](https://raw.githubusercontent.com/henkelmangroup/bader/main/examples/NaCl_wcore.tar.gz) crystal including core charges
  (vasp chgcar)
- [C<sub>2</sub>H<sub>4</sub>](https://raw.githubusercontent.com/henkelmangroup/bader/main/examples/C2H4a.tar.gz) molecule, orientation 1 (vasp
  chgcar)
- [C<sub>2</sub>H<sub>4</sub>](https://raw.githubusercontent.com/henkelmangroup/bader/main/examples/C2H4b.tar.gz) molecule, orientation 2 (vasp
  chgcar)
- [H<sub>2</sub>O](https://raw.githubusercontent.com/henkelmangroup/bader/main/examples/H2O.tar.gz) molecule (gaussian cube)

## Note for VASP users

One major issue with the charge density (CHGCAR) files from the VASP
code is that they only contain the valance charge density. The Bader
analysis assumes that charge density maxima are located at atomic
centers (or at pseudoatoms). Aggressive pseudopotentials remove charge
from atomic centers where it is both expensive to calculate and
irrelevant for the important bonding properties of atoms.

VASP contains a module (aedens) which allows for the core charge to be
written out from PAW calculations. By adding the LAECHG=.TRUE. to the
INCAR file, the core charge is written to AECCAR0 and the valance charge
to AECCAR2. These two charge density files can be summed using the
[chgsum.pl](https://github.com/henkelmangroup/vtstscripts) script:


      chgsum.pl AECCAR0 AECCAR2

The total charge will be written to CHGCAR_sum.

The bader analysis can then be done on this total charge density file:


      bader CHGCAR -ref CHGCAR_sum

One finally note is that you need a fine fft grid to accurately
reproduce the correct total core charge. It is essential to do a few
calculations, increasing NG(X,Y,Z)F until the total charge is correct.

## Note for CASTEP users

Aaron Hopkinson and Dr. Matt Probert from the University of York have
provided a [den2vasp.tar.gz](https://github.com/henkelmangroup/bader/releases/latest/download/den2vasp.tar.gz) utility to
convert from the CASTEP charge density to the VASP CHGCAR format so that
it can be read in by the Bader analysis program.

## Authors

This program was written by [Andri
Arnaldsson](mailto:andri@theochem.org), [Wenjie
Tang](mailto:tangwenjie@utexas.edu), [Sam
Chill](mailto:samchill@gmail.com), [Wenrui
Chai](mailto:chaiwenrui@gmail.com), [Ryan
Anselm](mailto:ryan.a.anselm@gmail.com%20), and [Graeme
Henkelman](mailto:henkelman@utexas.edu).

Improvements to the original algorithm were developed by Ed Sanville
(Loughborough University, UK).

With contributions from: Johannes Voss (DTU), Erik McNellis (FHI), and
Matthew Dyer (Liverpool)

Multipole code added by: Sebastien Lebegue, Angyan Janos, and Emmanuel
Aubert (Institut Jean Barriol, Nan cy-University)

## References

- W. Tang, E. Sanville, and G. Henkelman [A grid-based Bader analysis algorithm without lattice bias](https://theory.cm.utexas.edu/henkelman/code/bader/download/tang09_084204.pdf), *J. Phys.: Condens. Matter* **21**, 084204 (2009).

- E. Sanville, S. D. Kenny, R. Smith, and G. Henkelman [An improved grid-based algorithm for Bader charge allocation](https://theory.cm.utexas.edu/henkelman/code/bader/download/sanville07_899.pdf), *J. Comp. Chem.* **28**, 899-908 (2007).
  
- G. Henkelman, A. Arnaldsson, and H. Jónsson, [A fast and robust algorithm for Bader decomposition of charge density](https://theory.cm.utexas.edu/henkelman/code/bader/download/henkelman06_354.pdf), *Comput. Mater. Sci.* **36**, 354-360 (2006).

- M. Yu and D. R. Trinkle, [Accurate and efficient algorithm for Bader charge integration](https://theory.cm.utexas.edu/henkelman/code/bader/download/yu11_064111.pdf), *J. Chem. Phys.* **134**, 064111 (2011).
