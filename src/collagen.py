import re

def period(seq):
    pattern = r'^(?:G[A-Z]{2})+|(?:G[A-Z]{2})+G|(?:g[A-Z]{2})+G[A-Z]$'
    if re.fullmatch(pattern, seq):
        return 0
    pattern =r'^(?:[A-Z]G[A-Z]{2})+|(?:[A-Z]G[A-Z]{2})+G|(?:[A-Z]G[A-Z]{2})+G[a-z]$'
    if re.fullmatch(pattern, seq):
        return 1
    pattern =r'^(?:[A-Z]{2}G[A-Z]{2})+|(?:[A-Z]{2}G[A-Z]{2})+g|(?:[A-Z]{2}G[A-Z]{2})+G[A-Z]$'
    if re.fullmatch(pattern, seq):
        return 2
    return -1
    
def P_pattern(seq):
    p=period(seq)
    if p<0:
        return (-1, -1)
    weak_P=0
    strong_P=0
    for i in range((p+2)%3, len(seq), 3):
        if seq[i] == 'P':
            strong_P+= 1
    for i in range((p+1)%3, len(seq), 3):
        if seq[i] == 'P':
            weak_P+= 1
    return (strong_P, weak_P)
    
def Pmask_distance(seq1, seq2):
    if len(seq1)!=len(seq2):
        return -1
    r=0
    for c1, c2 in zip(seq1, seq2):
        if c1=='P' and c2=='P':
            None
        elif c1=='P'or c2=='P':
            r=r+1
    return r
    
