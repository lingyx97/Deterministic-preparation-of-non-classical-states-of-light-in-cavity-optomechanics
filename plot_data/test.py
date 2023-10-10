import qutip as qtp
from pathlib import Path


data=qtp.qload(Path(__file__).parent/"plot678samples8")
print([d[-5].ptrace(0).diag()[2] for d in data])
print([d[1] for d in data])
if [0,0,0]:
    print(0)
# print([d[-5].ptrace(0).diag() for d in data])