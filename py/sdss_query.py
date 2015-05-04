#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SYNOPSIS

    TODO helloworld [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    TODO: Show some examples of how to use this script.

EXIT STATUS

    TODO: List exit codes

AUTHOR

    TODO: Name <name@example.org>

LICENSE

    This script is in the public domain, free from copyrights or restrictions.

VERSION

    $
"""

from __future__ import division, print_function

__all__ = ["Module Name"]
__version__ = "0.0.0"
__author__ = "Jonatan Selsing (jselsing@dark-cosmology.dk)"
__copyright__ = "Copyright 2014 Jonatan Selsing"


#
#

def test():
    """ Testing Docstring"""
    pass


if __name__ == '__main__':

    from astroquery.sdss import SDSS
    SDSS.get_images()
    # import json
    # import urllib2
    #
    # query_terms = dict()
    # # query_terms["ra"] = str(ob[0].header['RA'])+'d' #"185.1d"
    # # query_terms["dec"] = str(ob[0].header['DEC'])  #"56.78"
    # # query_terms["radius"] = "10.0"
    #
    # http://skyserver.sdss.org/dr12/en/tools/search/x_sql.aspx?
    # format=xml&cmd=select top 10 * from galaxy
    # url = "http://api.sdss3.org/spectrumQuery?" + '&'.join(["{0}={1}".format(key, value) for key, value in query_terms.items()])
    # print(url)
    # # make call to API
    # response = urllib2.urlopen(url)
    #
    # # read response, converting JSON format to Python list
    # matching_ids = json.loads(response.read())
    # # print(json.dumps(matching_ids, indent=4))
    #
    # # get the first id
    # spec_id = matching_ids[0]
    #
    # url = "http://api.sdss3.org/spectrum?id={0}&format=json".format(spec_id)
    #
    # response = urllib2.urlopen(url)
    # result = json.loads(response.read())
    # SDSS_spectrum = result[spec_id]
    #
    # wl_sdss = np.array(SDSS_spectrum["wavelengths"])
    # flux_sdss =  np.array(SDSS_spectrum["flux"])
    # z_sdss = (np.array(SDSS_spectrum["z"]))