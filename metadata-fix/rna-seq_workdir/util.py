import os

entry_dir = 'rna-seq_gnos_entry'
id_fixes_dir = 'rna-seq_id_fixes'

def write_file(flist, fn):
    file_dir = os.path.dirname(fn)
    if not os.path.exists(file_dir): os.makedirs(file_dir)
    with open(fn, 'w') as f:
        header = True  
        for r in flist:
            if header:
                f.write('\t'.join(r.keys()) + '\n')
                header = False 
            # make the list of output from dict
            line = []
            for p in r.keys():
                if isinstance(r.get(p), list):
                    line.append('|'.join(r.get(p)))
                elif isinstance(r.get(p), set):
                    line.append('|'.join(list(r.get(p))))
                elif r.get(p) is None:
                    line.append('')
                else:
                    line.append(str(r.get(p)))
            f.write('\t'.join(line) + '\n') 
