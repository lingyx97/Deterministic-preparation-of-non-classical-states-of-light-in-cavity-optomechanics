import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent))


for i in range(2,6):
    with open(str(Path(__file__).parent/('evolution{}.txt'.format(str(i)))), 'r') as file:
        data = file.read().replace('\n', ' ').replace("[","(").replace("]",")")
        data=data.replace("zz075NonCommutativeTimes","Nct")
        data=data.replace("zz080HermitianConjugate(a)","ad").replace("zz080HermitianConjugate(b)","bd")

    with open(str(Path(__file__).parent/('evolution{}PY.txt'.format(str(i)))), 'w') as file:
        file.write(data)
    

with open(str(Path(__file__).parent/'test.txt'), 'r') as file:
    data = file.read().replace('\n', ' ').replace("[","(").replace("]",")")
    data=data.replace("zz075NonCommutativeTimes","Nct")
    data=data.replace("zz080HermitianConjugate(a)","ad").replace("zz080HermitianConjugate(b)","bd")

with open(str(Path(__file__).parent/'test.txt'), 'w') as file:
    file.write(data)
