import sys

genetic_code={'UUU':'F','UUC':'F','UUA':'L','UUG':'L','CUU':'L','CUC':'L','CUA':'L','CUG':'L','AUU':'I','AUC':'I','AUA':'I','AUG':'M','GUU':'V','GUC':'V','GUA':'V','GUG':'V','UCU':'S','UCC':'S','UCA':'S','UCG':'S','CCU':'P','CCC':'P','CCA':'P','CCG':'P','ACU':'T','ACC':'T','ACA':'T','ACG':'T','GCU':'A','GCC':'A','GCA':'A','GCG':'A','UAU':'Y','UAC':'Y','UAA':'Stop','UAG':'Stop','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q','AAU':'N','AAC':'N','AAA':'K','AAG':'K','GAU':'D','GAC':'D','GAA':'E','GAG':'E','UGU':'C','UGC':'C','UGA':'Stop','UGG':'W','CGU':'R','CGC':'R','CGA':'R','CGG':'R','AGU':'S','AGC':'S','AGA':'R','AGG':'R','GGU':'G','GGC':'G','GGA':'G','GGG':'G'}

def translate(rna):
    rna=rna.upper().replace('T','U')
    p=[]
    for i in range(0,len(rna),3):
        c=rna[i:i+3]
        if len(c)<3: break
        aa=genetic_code.get(c,'')
        if aa=='Stop': break
        p.append(aa)
    return ''.join(p)

seq = sys.argv[1] if len(sys.argv)>1 else input("Enter RNA/DNA sequence: ")
print(translate(seq))
