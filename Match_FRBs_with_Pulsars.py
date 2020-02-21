import numpy as np
import astropy.coordinates as coord
import astropy.units as u
import argparse


class frbInformation():
    '''
    Gets information from the FRB cat.
    '''
    def __init__(self, filename):
        '''
        Reads the FRB cat from FRBCat.org and get info.
        
        Args:
        filename (str): the filename (path included) to the
                        frbcat file
        '''
        
        FRB_header = np.genfromtxt(filename,
                                   dtype=str,
                                   delimiter='","')[0]
        FRB_header_string = []
        for val in FRB_header:
            FRB_header_string.append(val)
        FRB_header_string = np.array(FRB_header_string)
        
        self.telescope_index = np.where(FRB_header_string=='telescope')[0][0]
        self.gl_index = np.where(FRB_header_string=='rop_gl')[0][0]
        self.gb_index = np.where(FRB_header_string=='rop_gb')[0][0]
        self.dm_index = np.where(FRB_header_string=='rmp_dm')[0][0]
        self.width_index = np.where(FRB_header_string=='rmp_width')[0][0]
        self.name_index = np.where(FRB_header_string=='\ufeff"frb_name')[0][0]
        
        self.FRB_info = np.genfromtxt(filename,
                                      skip_header=1,
                                      dtype=str,
                                      delimiter='","')
        
        dm = self.FRB_info[:, self.dm_index]
        DM = []
        for i, d in enumerate(dm):
            if '&' in d:
                val = d[:d.index('&')]
                val = float(val)
            else:
                val = float(d)
            DM.append(val)
        self.DM = np.array(DM)
        self.logDM = np.log10(DM)
        
        self.gl = self.FRB_info[:, self.gl_index].astype(float)
        self.gb = self.FRB_info[:, self.gb_index].astype(float)
        
        self.name = self.FRB_info[:, self.name_index].astype(str)
        
        self.gl_corrected = np.copy(self.gl)
        self.gl_corrected[np.where(self.gl_corrected>180)] = (self.gl_corrected[np.where(self.gl_corrected>180)] -
                                                              360)
        self.gl_rad = np.deg2rad(self.gl_corrected)
        self.gb_rad = np.deg2rad(self.gb)
        
        self.telescopes = self.FRB_info[:, self.telescope_index]
        self.unique_telescopes = np.unique(self.telescopes)
        
        self.telescope_colours = [u'#7cb6ff',
                                  u'#558B2F',
                                  'White',
                                  u'#8dff76',
                                  '#3F51B5',
                                  'Grey',
                                  '#F7567C',
                                  '#EAC5D8',
                                  'Black',
                                  '#FFE74C',
                                  '#854D27']


def get_beam_dimensions(DEC, beam_radius, hour_angle=0.0):
    '''
    Determines the height and width of the MeerKAT primary
    beam at a specific declination
    
    Args:
    
    DEC (float): the declination of the phase centre of the
                 MeerKAT pointing in decimal degrees.
    beam_radius (float): the radius of the beam in decimal
                         degrees if the telescope
                         was pointing directly up. i.e. the radius
                         of the beam if the beam were circular.
    *kwargs:
    hour_angle (float): the hour angle in decimal degrees. In the
                        case of the incoherent beam this value is
                        not important. Default: 0.0
    '''
    # Don't forget deg2rad and rad2deg!!!
    
    latitude = np.deg2rad(-30.721) # latitude of MeerKAT
    straight_up = np.where(np.deg2rad(DEC)==latitude)[0]

    # For a nice explanation of these equations
    # see: http://www.stargazing.net/kepler/altaz.html
    sin_boresight_alt = ((np.sin(np.deg2rad(DEC)) *
                          np.sin(latitude)) +
                         (np.cos(np.deg2rad(DEC)) *
                          np.cos(latitude) *
                          np.cos(np.deg2rad(hour_angle))))
    boresight_alt = np.arcsin(sin_boresight_alt)

    cos_boresight_azim = ((np.sin(np.deg2rad(DEC)) -
                           np.sin(boresight_alt) *
                           np.sin(latitude)) /
                          (np.cos(boresight_alt) *
                           np.cos(latitude)))
    cos_boresight_azim[cos_boresight_azim > 1] = 1
    cos_boresight_azim[cos_boresight_azim < -1] = -1

    # Don't forget to adjust the azimuth value
    # depending on the hour angle
    boresight_azim = np.arccos(cos_boresight_azim)
    if hour_angle >= 0.:
        boresight_azim = 2. * np.pi - boresight_azim
    
    A0 = (boresight_azim + np.deg2rad(beam_radius))
    sinRA0 = (-np.sin(A0)*np.sin(boresight_alt)/
              np.cos(np.deg2rad(DEC)))
    sinRA0[sinRA0 > 1] = 1
    sinRA0[sinRA0 < -1] = -1
    RA0 = np.arcsin(sinRA0)

    WIDTH = np.rad2deg(RA0)
    # If the telescope is pointing exactly
    # up the beam is circular
    WIDTH[straight_up] = beam_radius
    # The height of the beam is always the
    # same
    HEIGHT = np.ones(len(DEC))*beam_radius

    return np.abs(WIDTH), HEIGHT


def get_angle(ra0, dec0, ra1, dec1):
    '''
    Get the angle between the horizontal of an ellipse
    and the vertical (zero if horizontal)
    
    Args:
    ra0 (array): the ra in decimal degrees of the first
                 source
    dec0 (array): the dec in decimal degrees of the first
                  source
    ra1 (array): the ra in decimal degrees of the second
                 source
    dec1 (array): the dec in decimal degrees of the second
                  source
                  
    Returns:
    An array of angles (in radian) between the two sources.
    '''
    x = np.abs(ra0 - ra1)
    y = np.abs(dec0 - dec1)
    
    return np.arctan(y / x)


def get_separation(ra0, dec0, ra1, dec1):
    '''
    Gets the separation between two sources.
    
    Uses basic pythagoras to get the
    distance between two sources.
    
    Args:
    ra0 (array): the ra in decimal degrees of the first
                 source
    dec0 (array): the dec in decimal degrees of the first
                  source
    ra1 (array): the ra in decimal degrees of the second
                 source
    dec1 (array): the dec in decimal degrees of the second
                  source
                  
    Returns:
    An array of the separations (in degrees) between
    the sources.
    '''
    x = (ra0 - ra1)
    y = (dec0 - dec1)
    
    return np.sqrt(x**2 + y**2)


def get_ellipse_radius(a, b, angles):
    '''
    The radius of an ellipse at a specific angle.
    
    Args:
    a (array): the semi-major axis of the ellipse
    b (array): the semi-minor axis of the ellipse
    angles (array): the angles in radian
    
    Returns:
    An array of the radii at those angles with those
    axis radii
    '''
    top = a * b
    bottom = np.sqrt((a**2 * np.sin(angles)**2) + (b**2 * np.cos(angles)**2))
    
    return top / bottom


def get_line(frb_c,
             pulsar_c,
             pulsar_name,
             radius,
             ellipse_radius,
             band,
             h,
             w,
             frb_name,
             frb_telescope):
    '''
    Makes a list of the info needed to write to file.
    
    Args:
    frb_c (Astropy SkyCoord): the coordinates of the FRB
    pulsar_c (Astropy SkyCoord): the coordinates of the pulsar
    pulsar_name (str): the name of the pulsar
    radius (float): the distance in degrees between the FRB
                    and the pulsar
    ellipse_radius (float): the radius of the incoherent beam
                            at the angle between the FRB and pulsar
    band (str): the name of the band (e.g. l-band)
    h (float): the semi-minor axis of the incoherent beam
               in degrees
    w (float): the semi-major axis of the incoherent beam
               in degrees
    frb_name(str): the name of the FRB
    frb_telescope (str): the name of the telescope that 
                         detected the FRB
                         
    Returns:
    A list of the format:
            [frb_name,
            frb_telescope,
            frb_ra,
            frb_dec,
            frb_ra_deg,
            frb_dec_deg,
            pulsar_name,
            pulsar_ra,
            pulsar_dec,
            pulsar_ra_deg,
            pulsar_dec_deg,
            radius,
            ellipse_radius,
            band,
            h,
            w]
    '''

    frb_name = frb_name.strip('\""')
    frb_telescope = frb_telescope.capitalize()
    frb_ra = frb_c.fk5.ra.to_string(u.hour)
    frb_dec = frb_c.fk5.dec.to_string(u.degree,
                                      alwayssign=True)
    frb_ra_deg = frb_c.fk5.ra.deg
    frb_dec_deg = frb_c.fk5.dec.deg

    pulsar_ra = pulsar_c.fk5.ra.to_string(u.hour)
    pulsar_dec = pulsar_c.fk5.dec.to_string(u.degree,
                                            alwayssign=True)
    pulsar_ra_deg = pulsar_c.fk5.ra.deg
    pulsar_dec_deg = pulsar_c.fk5.dec.deg

    return [frb_name,
            frb_telescope,
            frb_ra,
            frb_dec,
            frb_ra_deg,
            frb_dec_deg,
            pulsar_name,
            pulsar_ra,
            pulsar_dec,
            pulsar_ra_deg,
            pulsar_dec_deg,
            radius,
            ellipse_radius,
            band,
            h,
            w]


def match_frbs_pulsars(pulsar_coords,
                       frb_coords,
                       frb_info,
                       pulsar_names,
                       l_beam_radius=1.11,
                       uhf_beam_radius=2.31,
                       circular=False):
    '''
    Find out whether an FRB is within the incoherent
    beam when observing a pulsar with MeerKAT.
    
    Works out whether any FRBs are inside the incoherent
    beam of MeerKAT at L-band and UHF-band when MeerKAT
    is observing a pulsar (i.e. the pulsar is at the phase
    centre). Save the informatio in an array.
    
    Args:
    pulsar_coords (array): An array of astropy SkyCoord objects
                           where each value gives the coordinates
                           of a pulsar in the ATNF pulsar catalogue
    frbs_coords (array):An array of astropy SkyCoord objects
                        where each value gives the coordinates
                        of an FRB from the FRBcat catalogue
    frb_info (object): the object created by frbInformation
                       by reading the FRBcat file
    pulsar_names (array): an array of pulsar names matching the
                          coordinates in pulsar_coords
    kwargs:
    l_beam_radius (float): the beam radius in degrees for L-band.
                           Default: 1.11
    uhf_beam_radius (float): the beam radius in degrees for
                             UHF-band. Default:2.31
    circular (bool): whether or not to assume the beam is
                     always ideally circular, or to calculate
                     the elliptical projection on the sky.
                     Default: False (i.e. assume elliptical)
    Returns:
    An array of information about the FRBs that are near
    enough to pulsars to be inside the incoherent beam.
    '''
    if circular:
        print('Assuming the beam is always circular')
    else:
        print('Calculating elliptical beams')
    
    # Get the ellipse dimensions for all ATNF objects
    # for both bands
    l_w, l_h = get_beam_dimensions(pulsar_coords.fk5.dec.deg,
                                   l_beam_radius)
    uhf_w, uhf_h = get_beam_dimensions(pulsar_coords.fk5.dec.deg,
                                       uhf_beam_radius)
    match_info = []
    for f, frb in enumerate(frb_coords):
        # FRB RA and DEC
        f_ra = frb.fk5.ra.deg
        f_dec = frb.fk5.dec.deg
        # Pulsar RAs and DECs
        p_ra = pulsar_coords.fk5.ra.deg
        p_dec = pulsar_coords.fk5.dec.deg

        # Standard info for this FRB
        frb_c = frb
        frb_name = frb_info.name[f]
        frb_telescope = frb_info.telescopes[f]

        # Work out whether the FRB is inside
        # the incoherent beam for any pulsars
        radii = get_separation(f_ra, f_dec, p_ra, p_dec)
        angles = get_angle(f_ra, f_dec, p_ra, p_dec)

        # L-band
        if circular:
            ellipse_radii = np.ones(len(pulsar_coords))*l_beam_radius
        else:
            ellipse_radii = get_ellipse_radius(l_w, l_h, angles)
        inside = np.where(ellipse_radii - radii >= 0)[0]
        # Find FRBs and make a row of info about them
        if len(inside) > 0:
            band = 'l-band'
            for val in inside:
                pulsar_c = pulsar_coords[val]
                pulsar_name = pulsar_names[val]
                radius = radii[val]
                ellipse_radius = ellipse_radii[val]
                h = l_h[val]
                w = l_w[val]

                match_info.append(get_line(frb_c,
                                           pulsar_c,
                                           pulsar_name,
                                           radius,
                                           ellipse_radius,
                                           band,
                                           h,
                                           w,
                                           frb_name,
                                           frb_telescope))
        # UHF-band
        if circular:
            ellipse_radii = np.ones(len(pulsar_coords))*uhf_beam_radius
        else:
            ellipse_radii = get_ellipse_radius(uhf_w, uhf_h, angles)
        inside = np.where(ellipse_radii - radii >= 0)[0]
        # Find FRBs and make a row of info about them
        if len(inside) > 0:
            band = 'uhf-band'
            for val in inside:
                pulsar_c = pulsar_coords[val]
                pulsar_name = pulsar_names[val]
                radius = radii[val]
                ellipse_radius = ellipse_radii[val]
                h = uhf_h[val]
                w = uhf_w[val]

                match_info.append(get_line(frb_c,
                                           pulsar_c,
                                           pulsar_name,
                                           radius,
                                           ellipse_radius,
                                           band,
                                           h,
                                           w,
                                           frb_name,
                                           frb_telescope))
    return np.array(match_info, dtype=str)


def get_pulsar_info(filename, deg_cut=44.):
    '''
    Reads the ATNF catalogue file.
    
    Args:
    filename (str): the name of the csv of
                    the ATNF pulsar catalogue
    kwargs:
    deg_cut (float): the maximum declination
                     in degrees that you're telescope can
                     observe. Default: 44
    Returns:
    An array of pulsar coordinates as SkyCoord
    objects, and an array of pulsar names
    (matching the pulsar coordinates).
    '''
    # Read in the header to get the
    # columns. This is because you'll
    # probably want to update this one over
    # time and the columns you get
    # might change. It's safer to just
    # download all the columns that are
    # available and select for the ones
    # you actually want
    with open(filename, 'r') as f:
        hdr = f.readline().strip()
    hdr = hdr.split(';')
    hdr = np.array(hdr)
    # Use the header to find the index of
    # the columns that you're actually
    # interested in
    jname_index = np.where(hdr=='PSRJ')[0][0]
    bname_index = np.where(hdr=='NAME')[0][0]
    ra_index = np.where(hdr=='RAJ')[0][0]
    dec_index = np.where(hdr=='DECJ')[0][0]

    # This catalogue can be read using
    # genfromtxt, it has a couple lines
    # for header, and is semi-colon
    # separated
    values = np.genfromtxt(filename,
                           delimiter=';',
                           skip_header=2,
                           dtype=str)
    # Get the columns and info
    # that you're actuall interested
    # in using the column indices found
    # by looking at the header
    ras_hms = values[:, ra_index]
    decs_dms = values[:, dec_index]
    pulsar_coords = coord.SkyCoord(ra=ras_hms, dec=decs_dms,
                                   frame='icrs',
                                   unit=(u.hourangle, u.deg))
    # You'll want the full pulsar names,
    # including "PSR" and if it's a
    # B or J name (or both) for when
    # you look the objects up in SIMBAD
    bnames = values[:, bname_index]
    bnames = np.char.replace(bnames, 'J', 'PSR J')
    bnames = np.char.replace(bnames, 'B', 'PSR B')
    pulsar_names = np.copy(bnames)

    # Only pulsars below 44 deg declination
    pulsar_decs = np.where(pulsar_coords.fk5.dec.deg<deg_cut)[0]

    pulsar_coords = pulsar_coords[pulsar_decs]
    pulsar_names = pulsar_names[pulsar_decs]
    
    return pulsar_coords, pulsar_names


def write_match_information(pulsar_file,
                            frb_file,
                            output_file,
                            l_beam_radius=1.11,
                            uhf_beam_radius=2.312,
                            deg_cut=44.0,
                            circular=False):
    '''
    Find FRB pulsar matches and write the information
    to file.
    
    Args:
    pulsar_file (str): the name and path of the ATNF
                       pulsar catalogue file
    frb_file (str): the name and path of the FRBcat file
    output_file (str): the name and path of the output file
    
    kwargs:
    l_beam_radius (float): the beam radius in degrees for L-band.
                           Default: 1.11
    uhf_beam_radius (float): the beam radius in degrees for
                             UHF-band. Default:2.31
    circular (bool): whether or not to assume the beam is
                     always ideally circular, or to calculate
                     the elliptical projection on the sky.
    deg_cut (float): the maximum declination
                     in degrees that you're telescope can
                     observe. Default: 44
    Returns:
    An array of information, the same information
    that is written to output_file.
    '''
    frb_info = frbInformation(frb_file)

    frb_coords = coord.SkyCoord(frb_info.gl,
                                frb_info.gb,
                                frame='galactic',
                                unit='deg')
    pulsar_coords, pulsar_names = get_pulsar_info(pulsar_file,
                                                   deg_cut=deg_cut)
    
    match_info = match_frbs_pulsars(pulsar_coords,
                                    frb_coords,
                                    frb_info,
                                    pulsar_names,
                                    l_beam_radius=l_beam_radius,
                                    uhf_beam_radius=l_beam_radius,
                                    circular=circular)
    
    # Write it all to a file

    header = ('FRB name, FRB telescope, '
              'FRB RA, FRB DEC, '
              'FRB RA(deg), FRB DEC (deg), '
              'Pulsar name, '
              'Pulsar RA, Pulsar DEC, '
              'Pulsar RA (deg), Pulsar DEC (deg), '
              'Separation (deg), Ellipse radius (deg), '
              'Band, '
              'Beam height, Beam width')
    np.savetxt(output_file,
               match_info, header=header,
               delimiter=',', fmt='%s')

    print('---DONE---')
    print('{} written to file'.format(output_file))
    
    return match_info


if __name__ in '__main__':
    parser = argparse.ArgumentParser(description=('Find whether any FRBs are in the '
                                                  'MeerKAT FoV of any pulsars'))

    parser.add_argument('-F', action='store',
                        dest='FRB_file',
                        help=('The name and path of the FRB catalogue file '
                              '(file downloaded from FRBcat.org)'))
    parser.add_argument('-P',action='store',
                        dest='PSR_file',
                        help=('The name and path of the ATNF pulsar catalogue '
                              'file (file downloaded from '
                              'https://www.atnf.csiro.au/research/pulsar/psrcat/)'))
    parser.add_argument('--o', action='store',
                        dest='output_file',
                        help=('The name and path of the final csv file. '
                              'Default: Pulsar_FRB_Matches.csv'),
                        default='Pulsar_FRB_Matches.csv')

    parser.add_argument('--l_beam', action='store',
                        dest='l_beam',
                        help=('The radius of the L-band beam in degrees. '
                              'Default: 1.11'), default=1.11)
    parser.add_argument('--uhf_beam', action='store',
                        dest='uhf_beam',
                        help=('The radius of the UHF-band beam in degrees. '
                              'Default: 2.31'),
                        default=2.31)
    parser.add_argument('--max_dec', action='store',
                        dest='dec_cut',
                        help=('The maximum declination (in degrees) you can observe, '
                              'excludes any pulsars/FRBs above that declination. '
                              'Default: 44.0'),
                        default=44.0)
    parser.add_argument('--circular', action='store_true',
                        default=False, dest='circular',
                        help=('Use the beam shape as the theoretical circular '
                              'beam no matter where the telescope is pointing. '
                              'Default: False (i.e. assume elliptical beam shape, '
                              'not circular).'))


    results = parser.parse_args()


    write_info = write_match_information(results.PSR_file,
                                         results.FRB_file,
                                         results.output_file,
                                         l_beam_radius=float(results.l_beam),
                                         uhf_beam_radius=float(results.uhf_beam),
                                         deg_cut=float(results.dec_cut),
                                         circular=results.circular)
