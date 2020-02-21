# FRBs_in_Pulsar_Obs

### Find FRBs that are within the primary beam of MeerKAT when a pulsar is at the phase centre

Pulsars are regularly observed by MeerKAT teams such as MeerTIME and TRAPUM. The MeerTRAP team
commensally searches for Fast Radio Bursts (FRBs) with MeerKAT projects. This means that MeerTRAP
observes the locations of some FRBs when piggybacking on MeerTIME and TRAPUM projects.

The iPython notebook ad python script in this repo use the positions of pulsars from the
<a href="https://www.atnf.csiro.au/research/pulsar/psrcat/">ATNF</a> catalogue to calculate the
size of the MeerTRAP beam. This takes into account the ellipticity of the beam at different
declinations. It then calculates the separation between all of the pulsars and each FRB in the
<a href="http://frbcat.org/">FRBcat</a> and decides whether the FRB is inside the MeerKAT beam when
each pulsar is at the phase centre. You can give the code the size of the MeerKAT beam assuming it's
circular, or it will assume the radius of the beam to be 1.11 degrees at L-band and 2.31 degrees for
the UHF band.

Finally, the code writes out a CSV file that has columns with:
<ul>
  <li>FRB name</li>
  <li>FRB telescope</li>
  <li>FRB RA</li>
  <li>FRB DEC</li>
  <li>FRB RA (deg)</li>
  <li>FRB DEC (deg)</li>
  <li>Pulsar name</li>
  <li>Pulsar RA</li>
  <li>Pulsar DEC</li>
  <li>Pulsar RA (deg)</li>
  <li>Pulsar DEC (deg)</li>
  <li>Separation (deg)</li>
  <li>Ellipse radius (deg)</li>
  <li>Band</li>
  <li>Beam height</li>
  <li>Beam width</li>
</ul>

Part of this code was originally devoloped by
<a href="https://github.com/BezuidenhoutMC/MosaicUtils">BezuidenhoutMC</a> as part of his Mosaic code.
I adapted <a href="https://github.com/BezuidenhoutMC/MosaicUtils/blob/master/Primary_plotter.py">Primary_plotter.py</a>
to calculate the shape of the primary beam on the sky.

## What's in this respository

<ul>
  <li>ATNF.csv</li>
  <li>frbcat_20200129.csv</li>
  <li>Match_FRBs_with_Pulsars.ipynb</li>
  <li>Match_FRBs_with_Pulsars.py/li>
</ul>

#### ATNF.csv
This is a download of the <a href="https://www.atnf.csiro.au/research/pulsar/psrcat/">ATNF</a> catalogue.

#### frbcat_20200129.csv
This is a download of <a href="http://frbcat.org/">FRBcat</a> on 2020 February 29 using the <tt>download</tt> button
on the website.

#### Match_FRBs_with_Pulsars.ipynb
The jupyter notebook for finding FRBs in the MeerKAT beam.

#### Match_FRBs_with_Pulsars.py
The python script for finding FRBs in the MeerKAT beam.

The usage of Match_FRBs_with_Pulsars.py is:
Match_FRBs_with_Pulsars.py<br> [-h] [-F FRB_FILE] [-P PSR_FILE]<br>
                                  [--o OUTPUT_FILE] [--l_beam L_BEAM]<br>
                                  [--uhf_beam UHF_BEAM] [--max_dec DEC_CUT]<br>
                                  [--circular]

Find whether any FRBs are in the MeerKAT FoV of any pulsars

optional arguments:<br>
  -h, --help<br>           show this help message and exit<br>
  -F FRB_FILE<br>          The name and path of the FRB catalogue file (file
                       downloaded from FRBcat.org)<br>
  -P PSR_FILE<br>          The name and path of the ATNF pulsar catalogue file
                       (file downloaded from
                       https://www.atnf.csiro.au/research/pulsar/psrcat/)<br>
  --o OUTPUT_FILE<br>      The name and path of the final csv file. Default:
                       Pulsar_FRB_Matches.csv<br>
  --l_beam L_BEAM<br>      The radius of the L-band beam in degrees. Default: 1.11<br>
  --uhf_beam UHF_BEAM<br>  The radius of the UHF-band beam in degrees. Default:
                       2.31<br>
  --max_dec DEC_CUT<br>    The maximum declination (in degrees) you can observe,
                       excludes any pulsars/FRBs above that declination.
                       Default: 44.0<br>
  --circular<br>           Use the beam shape as the theoretical circular beam no
                       matter where the telescope is pointing. Default: False
                       (i.e. assume elliptical beam shape, not circular).<br>
                       
Please use <tt>python Match_FRBs_with_Pulsars.py -h</tt> to print the help message.
                       
For example:
<tt>python Match_FRBs_with_Pulsars.py -F frbcat_20200129.csv -P ATNF.csv --o Matching_Pulsars_and_FRBs.csv</tt>

<tt>--circular</tt> sets the beams to circular beams, i.e. does not take into
account the stretching of the beam at different declinations.
