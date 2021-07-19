from collections import defaultdict
from itertools import chain
from operator import methodcaller
import pickle as pkl
from tqdm import tqdm

all_pdbs = open('pdb_lists/unique_pdbs.txt').read().split(',')[:3]
all_coords = [pkl.load(open('data/extract_out/'+pdb+'.pkl', 'rb')) for pdb in all_pdbs]

combined_data = defaultdict(list)

dict_items = map(methodcaller('items'), all_coords)
for k, v in tqdm(chain.from_iterable(dict_items)):
        combined_data[k].extend(v)

pkl.dump(combined_data, open('data/raw_contacts/all_contacts.pkl', 'wb'))

