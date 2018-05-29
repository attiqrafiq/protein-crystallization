import re
import math
import numpy as np
from sklearn.svm import SVC
from sklearn import metrics
from sklearn.model_selection import cross_val_predict

# extract feature GE
def GE(protein_seq,partition_n):
    re_seq1 = re.compile('A|V|L|I|M|C')
    protein_seq = re_seq1.sub('1',protein_seq)
    re_seq2 = re.compile('F|W|Y|H')
    protein_seq = re_seq2.sub('2', protein_seq)
    re_seq3 =  re.compile('S|T|N|Q')
    protein_seq = re_seq3.sub('3', protein_seq)
    re_seq4 = re.compile('K|R')
    protein_seq = re_seq4.sub('4', protein_seq)
    re_seq5 = re.compile('D|E')
    protein_seq = re_seq5.sub('5', protein_seq)
    re_seq6 = re.compile('G|P')
    protein_seq = re_seq6.sub('6', protein_seq)

    Len = len(protein_seq)

    partitions = ('123', '124', '125', '126', '134', '135', '136', '145', '146', '156')
    binary_seq = []
    for i in range(10):
        CS_pt = ''
        pt = partitions[i]
        for j in range(Len):
            AA_str = protein_seq[j]
            if pt.find(AA_str) != -1:
                CS_pt += '1'
            else:
                CS_pt += '0'
        binary_seq.append(CS_pt)

    Feature = []
    for i in range(10):
        ss = binary_seq[i]
        for j in range(partition_n):
            l = math.floor(Len/float(partition_n))
            ll = (j+1) * l
            if (l*(j+2)) > Len:
                substr = ss[:]
            else:
                substr = ss[0:int(ll)]
            num_0 = substr.count('0')
            frequence_0 = float(num_0) / len(substr)
            Feature.append(frequence_0)
            num_1 = substr.count('1')
            frequence_1 = float(num_1) / len(substr)
            Feature.append(frequence_1)
            num_01 = substr.count('01')
            num_10 = substr.count('10')
            frequence_0_1 = float(num_01 + num_10) / len(substr)
            Feature.append(frequence_0_1)
    return Feature
# extract feature PsePSSM
def PsePSSM(P, lg):
    n = P.shape[0]
    if n < lg:
        row = [[0 for i in range(20)] for j in range(lg - n + 1)]
        row = np.array(row)
        P = np.row_stack((P, row))
        n = lg + 1

    P = np.array(P, dtype=float)
    ME = np.mean(P, 1)
    st = np.std(P, 1)
    for i in range(len(st)):
        if st[i] == 0.0:
            st[i] = 1.0
    for i in range(n):
        P[i] = (P[i]-ME[i])/st[i]
    feature = []
    for i in range(20):
        feature.append(sum(P[:, i]) / n)
    for lag in range(1, lg+1):
        a = 0.0
        for j in range(20):
            for k in range(n-lag):
                a = a+pow((P[k][j]-P[k+lag][j]), 2)
            a = a/(n-lag)
            feature.append(a)
    return feature
# extract feature AVBlock_PSSM and AVBlock_SS
def AVBlock_PSSM_SS(matrix):
    feature = []
    n = matrix.shape[0]
    m = matrix.shape[1]
    print 'm', m
    size_2 = int(math.floor(n/2.0))
    size_3 = int(math.floor(n/3.0))
    size_4 = int(math.floor(n/4.0))
    size_5 = int(math.floor(n/5.0))
    #divide 1
    for i in range(m):
        feature.append(float(sum(matrix[:, i])) /  n)
    print len(feature)
    #divide 2
    for i in range(2):
        if i != 1:
            block = matrix[size_2*i:size_2*(i+1), :]
            for j in range(m):
                feature.append(float(sum(block[:, j])) / size_2)
        else:
            block = matrix[size_2*i:n, :]
            for j in range(m):
                feature.append(float(sum(block[:, j])) / (n-size_2*1))
    print len(feature)
    #divide 3
    for i in range(3):
        if i != 2:
            block = matrix[size_3*i:size_3*(i+1), :]
            for j in range(m):
                feature.append(float(sum(block[:, j])) / size_3)
        else:
            block = matrix[size_3*i:n, :]
            for j in range(m):
                feature.append(float(sum(block[:, j])) / (n-size_3*2))
    print len(feature)
    #divide 4
    for i in range(4):
        if i != 3:
            block = matrix[size_4*i:size_4*(i+1), :]
            for j in range(m):
                feature.append(float(sum(block[:, j])) / size_4)
        else:
            block = matrix[size_4*i:n, :]
            for j in range(m):
                feature.append(float(sum(block[:, j])) / (n-size_4*3))
    print len(feature)
    #divide 5
    for i in range(5):
        if i != 4:
            block = matrix[size_5*i:size_5*(i+1), :]
            for j in range(m):
                feature.append(float(sum(block[:, j])) / size_5)
        else:
            block = matrix[size_5*i:n, :]
            for j in range(m):
                feature.append(float(sum(block[:, j])) / (n-size_5*4))
    print len(feature)
    return feature
# extract feature Protscale
def pipei(protein, p, i):
    if protein == 'A':
        return p[i][0]
    if protein == 'R':
        return p[i][1]
    if protein == 'N':
        return p[i][2]
    if protein == 'D':
        return p[i][3]
    if protein == 'C':
        return p[i][4]
    if protein == 'Q':
        return p[i][5]
    if protein == 'E':
        return p[i][6]
    if protein == 'G':
        return p[i][7]
    if protein == 'H':
        return p[i][8]
    if protein == 'I':
        return p[i][9]
    if protein == 'L':
        return p[i][10]
    if protein == 'K':
        return p[i][11]
    if protein == 'M':
        return p[i][12]
    if protein == 'F':
        return p[i][13]
    if protein == 'P':
        return p[i][14]
    if protein == 'S':
        return p[i][15]
    if protein == 'T':
        return p[i][16]
    if protein == 'W':
        return p[i][17]
    if protein == 'Y':
        return p[i][18]
    if protein == 'V':
        return p[i][19]
    print protein
def protscal(protein_seq, p, lg):
    length = len(protein_seq)
    feature = []
    for i in range(57):
        num_seq = []
        for j in range(length):
            if j < 3 or j > length-4:
                num_seq.append(pipei(protein_seq[j],p,i))
            else:
                num_seq.append((pipei(protein_seq[j-3],p,i)*0.1+pipei(protein_seq[j-2],p,i)*0.4+pipei(protein_seq[j-1],p,i)*0.7+pipei(protein_seq[j],p,i)+
                               pipei(protein_seq[j+1],p,i)*0.7+pipei(protein_seq[j+2],p,i)*0.4+pipei(protein_seq[j+3],p,i)*0.1)/3.4)
        feature1 = []
        for lag in range(1, lg+1):
            a = 0.0
            for k in range(length-lag):
                a = a+(num_seq[k]*num_seq[k+lag])
            a = a/(length-lag)
            feature1.append(a)
        feature.extend(feature1)
    return feature


#five-fold
x = trainx
y = trainy
print x.shape, y.shape
model=SVC(C=2,gamma=0.0156)
predict_label=cross_val_predict(model,x,y,cv=5)
tn, fp, fn, tp = metrics.confusion_matrix(y, predict_label).ravel()
print tn, fp, fn, tp
print 'ACC:', metrics.accuracy_score(y, predict_label)
print 'MCC:', metrics.matthews_corrcoef(y, predict_label)
print 'special:', float(tn) / (tn + fp)
print 'SN:', float(tp) / (tp + fn)

#independent test
x_test = testx
y_test = testy
print x_test.shape, y_test.shape
model=SVC(C=2,gamma=0.0156)
model.fit(x, y)
predict_label = model.predict(x_test)
print y_test.shape, predict_label.shape
tn, fp, fn, tp = metrics.confusion_matrix(y_test, predict_label).ravel()
print tn, fp, fn, tp
print 'ACC:', metrics.accuracy_score(y_test, predict_label)
print 'MCC:', metrics.matthews_corrcoef(y_test, predict_label)
print 'special:', float(tn) / (tn + fp)
print 'SN:', float(tp) / (tp + fn)
