#!/usr/bin/env python
# encoding: utf-8

'''
Original code by boyko from YeoLab/FLARE
Repository: [GitHub URL]
Modifications made by JeremyTse on Jan 10, 2025
Changes: Adapted for MAPIT-seq data, added functionality for reditools RNA-editing analysis.
'''

import sys
import os
import re
import pandas as pd
from scipy.special import betainc

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2016-07-13'
__updated__ = '2016-07-13'

BASE=['A','C','G','T']

COMPLEMENT = {
    'C': 'G',
    'A': 'T',
    'G': 'C',
    'T': 'A'
}

DEAMINATION = {
    'C': 'T',
    'A': 'G',
    'T': 'C',
    'G': 'A'
}

def rank_edits(alpha, beta, cov_margin, input_eff, outfile):
    df_redi_result=pd.read_csv(input_eff,sep="\t",header=None,names=["Region","Position","Reference","Strand","Coverage","MeanQ","BaseCount[A,C,G,T]","AllSubs","Frequency","gCoverage","gMeanQ","gBaseCount[A,C,G,T]","gAllSubs","gFrequency"])
    df_redi_result['START']=df_redi_result['Position']-1
    df_redi_result['ALT']=df_redi_result['Reference'].replace(DEAMINATION)
    df_redi_result['EDIT']=df_redi_result.apply(lambda x:int(x['BaseCount[A,C,G,T]'][1:-1].split(", ")[BASE.index(x["ALT"])]),axis=1)
    df_redi_result['EDITRate']=df_redi_result['EDIT']/df_redi_result['Coverage']
    strand = "+" if '.fwd.' in os.path.basename(input_eff) else "-"
    df_redi_result['STRAND']=strand
    df_redi_result['CONFIDENCE']=1-df_redi_result.apply(lambda x:betainc(x['EDIT'] + alpha, x['Coverage'] + beta, cov_margin),axis=1)
    df_redi_result['NAME']=df_redi_result["Coverage"].astype(str)+"|"+df_redi_result['Reference']+">"+df_redi_result['ALT']+"|"+df_redi_result['EDITRate'].astype(str)
    df_redi_result[["Region","START","Position","NAME","CONFIDENCE","STRAND"]].to_csv(outfile,sep="\t",header=None,index=False)

def main(argv=None):  # IGNORE:C0111
    """
    Command line options.
    :param argv:
    :return:
    """

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (
        program_version, program_build_date
    )
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s
  Repository: [GitHub URL]
  Modifications made by JeremyTse on Jan 10, 2025
  Changes: Adapted for MAPIT-seq data, added functionality for reditools RNA-editing analysis.
  Created by user_name on %s.
  Copyright 2016 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    # Setup argument parser
    parser = ArgumentParser(
        description=program_license,
        formatter_class=RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=program_version_message
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_eff",
        help="variant reditools2 result file for filtering",
        required=True
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_bed",
        help="output noSNP (bed) file with known snps removed.",
        required=True
    )
    parser.add_argument(
        "-c", "--cov_margin",
        dest="cov_margin",
        help="minimum edit fraction",
        required=False,
        default=0.01,
        type=float
    )
    parser.add_argument(
        "-a", "--alpha",
        dest="alpha",
        help="alpha parameter",
        required=False,
        default=0,
        type=int
    )
    parser.add_argument(
        "-b", "--beta",
        dest="beta",
        help="beta parameter",
        required=False,
        default=0,
        type=int
    )

    # Process arguments
    args = parser.parse_args()
    input_eff = args.input_eff
    outfile = args.output_bed
    cov_margin = args.cov_margin
    alpha = args.alpha
    beta = args.beta

    rank_edits(alpha, beta, cov_margin, input_eff, outfile)
    return 0

if __name__ == "__main__":
    sys.exit(main())