import os
import glob
import sys
import pandas as pd

in_directory=sys.argv[1]
prefix=sys.argv[2]
outname = prefix + '.xlsx'

search= in_directory  + '/*tab'
print(search)
writer = pd.ExcelWriter(outname, engine='xlsxwriter')

for f in sorted(glob.glob(search)):
    print(f)
    df = pd.read_csv(f, sep="\t")
    df.to_excel(writer, sheet_name=os.path.basename(f)[:31])

writer.save()
