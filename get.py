''' author: samtenka
    change: 2022-05-07
    create: 2022-05-07
    descrp: 
    to use: Run `python3 get.py` 
            Enter a (space-separated) list of tags, e.g. 'algtop diftop'
            The interface should look like:
                AVAILABLE TAGS: ...
                SEARCH FOR: algtop
                0 -- bott.1995.differential-forms-topology.djvu
                1 -- bott.1995.forms-algebraic-topology.djvu
                2 -- hilton.1953.introduction-homotopy-theory.djvu
                3 -- milnor.1974.characteristic-classes.djvu
                4 -- sato.1999.algebraic-topology-intuitively.djvu
                5 -- switzer.1975.homotopy-homology.djvu
                6 -- ueno.2003.a-mathematical-gift.djvu
                7 -- vick.1994.homology-algebraic-topology.djvu
                8 -- hu.1955.axiomatic-homotopy.pdf
                open? 
            Finally, enter a number indicating which result to open, e.g.:
                7 -- vick.1994.homology-algebraic-topology.djvu
                8 -- hu.1955.axiomatic-homotopy.pdf
                open? 7 
            And it should open in a doc reader!
'''

with open('label-descriptions') as f:
    labels = [parts[0]
              for l in f.read().split('\n') if l
              for parts in [l.strip().split()]  ]
    labels = sorted(labels)

with open('annotations') as f:
    tags_by_filenm = {'../{}s/{}'.format(parts[0],parts[-1]):set(parts[1:-1])
                      for l in f.read().split('\n') if l
                      for parts in [l.split()]          }

def search_for(query):
    desired_tags = set(query.strip().split())
    return [fnm for fnm, tags in tags_by_filenm.items() if desired_tags.issubset(tags)] 

import os
if __name__=='__main__':
    print('AVAILABLE TAGS: {}'.format(', '.join(labels)))
    while True:
        query = input('SEARCH FOR: ')
        matches = search_for(query)
        for i,fnm in enumerate(matches):
            print('{:3d} -- {}'.format(i,fnm.split('/')[-1]))
        query = input('open? ')
        if query:
            #os.system('xdg-open {}'.format(matches[int(query)]))
            os.system('evince {}'.format(matches[int(query)]))
