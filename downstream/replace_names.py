#! /usr/bin/env python3

import sys
import re
from get_mapping import get_mapping

if __name__ == "__main__":

    org_mapping,chr_mapping = get_mapping()
    with open(sys.argv[1]) as f:
        txt = f.read()
    js = [x for x in list(reversed(sorted(chr_mapping,key=len))) if re.match(r'\d+org',x)]
    for j in js:
        ks = [x for x in list(reversed(sorted(chr_mapping[j],key=len))) if 'orgchr' in x]
        for k in ks:
            txt = re.sub(k,chr_mapping[j][k],txt)
    ks = [x for x in list(reversed(sorted(org_mapping,key=len))) if re.match(r'\d+org',x)]
    for k in ks:
        txt = re.sub(k,org_mapping[k],txt)
    with open(sys.argv[2],'w') as f:
        f.write(txt)
